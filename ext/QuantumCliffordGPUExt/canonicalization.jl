using CUDA
using LinearAlgebra
using QuantumClifford  
const MIN_BLOCKS = 36

function swap_columns_and_phases_kernel!(d_mat::CuDeviceArray{T,2}, 
                                        d_phases::CuDeviceArray{UInt8,1}, 
                                        col1::Integer, col2::Integer, N::Integer) where {T}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx <= N
        temp = d_mat[idx, col1]
        d_mat[idx, col1] = d_mat[idx, col2]
        d_mat[idx, col2] = temp
    end
    if idx == 1
        temp_phase = d_phases[col1]
        d_phases[col1] = d_phases[col2]
        d_phases[col2] = temp_phase
    end
    return
end

function gpu_swap_columns_and_phases!(d_mat::CuArray{T,2}, d_phases::CuVector{UInt8}, 
                                     col1::Integer, col2::Integer) where {T}
    N = size(d_mat, 1)
    threads = 256 
    blocks = max(cld(N, threads), MIN_BLOCKS)
    @cuda threads=threads blocks=blocks swap_columns_and_phases_kernel!(
        d_mat, d_phases, col1, col2, N)
    return
end

function swap_columns_kernel!(d_mat::CuDeviceArray{T,2}, 
                             col1::Integer, col2::Integer, N::Integer) where {T}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx <= N
        temp = d_mat[idx, col1]
        d_mat[idx, col1] = d_mat[idx, col2]
        d_mat[idx, col2] = temp
    end
    return
end

function gpu_swap_columns!(d_mat::CuArray{T,2}, col1::Integer, col2::Integer) where {T}
    N = size(d_mat, 1)
    threads = 256
    blocks = max(cld(N, threads), MIN_BLOCKS)
    @cuda threads=threads blocks=blocks swap_columns_kernel!(d_mat, col1, col2, N)
    return
end

@inline function gpu_compute_phase(d_mat::CuDeviceArray{T,2}, pivot_col::Integer, 
                                  target_col::Integer, nblocks::Integer)::UInt8 where {T}
    cnt1 = zero(T)
    cnt2 = zero(T)
    @inbounds for i in 1:nblocks
        x1 = d_mat[i, pivot_col]
        z1 = d_mat[nblocks + i, pivot_col]
        x2 = d_mat[i, target_col]
        z2 = d_mat[nblocks + i, target_col]
        newx = x1 ⊻ x2
        newz = z1 ⊻ z2
        x1z2 = x1 & z2
        anti = (x2 & z1) ⊻ x1z2
        inner = (cnt1 ⊻ newx) ⊻ newz ⊻ x1z2
        cnt2 = cnt2 ⊻ (inner & anti)
        cnt1 = cnt1 ⊻ anti
    end
    rcnt1 = count_ones(cnt1)
    rcnt2 = count_ones(cnt2)    
    return UInt8((rcnt1 ⊻ (rcnt2 << 1)) & 0x3)
end

function pivot_search_kernel!(d_mat::CuDeviceArray{T,2}, pivot_row::Integer, bit::T, 
                             current_col::Integer, r::Integer, d_result) where {T}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x + current_col - 1
    
    if idx <= r
        if (d_mat[pivot_row, idx] & bit) != zero(T)
            CUDA.atomic_min!(pointer(d_result), Int32(idx))
        end
    end
    return
end

function gpu_pivot_search!(d_mat::CuArray{T,2}, pivot_row::Integer, bit::T, 
                          current_col::Integer, r::Integer)::Int where {T}
    n_cols = r - current_col + 1
    if n_cols <= 0
        return 0
    end
    d_result = CUDA.fill(Int32(r + 1), 1)
    threads = min(256, n_cols)
    blocks = cld(n_cols, threads)
    @cuda threads=threads blocks=blocks pivot_search_kernel!(
        d_mat, pivot_row, bit, current_col, r, d_result)
    result = Array(d_result)[1]
    return result > r ? 0 : Int(result)
end

function update_kernel!(d_mat::CuDeviceArray{T,2}, d_phases::CuDeviceArray{UInt8,1},
                       current_col::Integer, pivot_row::Integer, bit::T,
                       r::Integer, total_blocks::Integer, nblocks::Integer, 
                       update_phase::Bool) where {T}
    s_pivot = @cuDynamicSharedMem(T, total_blocks)
    tid = threadIdx().x
    for idx in tid:blockDim().x:total_blocks
        s_pivot[idx] = d_mat[idx, current_col]
    end
    sync_threads()
    m = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if m <= r && m != current_col
        if (d_mat[pivot_row, m] & bit) != zero(T)
            cnt1 = zero(T)
            cnt2 = zero(T)
            @inbounds for i in 1:nblocks
                x1 = s_pivot[i]
                z1 = s_pivot[nblocks + i]
                x2 = d_mat[i, m]
                z2 = d_mat[nblocks + i, m]
                newx = x1 ⊻ x2
                newz = z1 ⊻ z2
                x1z2 = x1 & z2
                anti = (x2 & z1) ⊻ x1z2
                inner = (cnt1 ⊻ newx) ⊻ newz ⊻ x1z2
                cnt2 = cnt2 ⊻ (inner & anti)
                cnt1 = cnt1 ⊻ anti
            end
            rcnt1 = count_ones(cnt1)
            rcnt2 = count_ones(cnt2)
            rcnt2_shifted = rcnt2 << 1
            extra_phase_val = rcnt1 ⊻ rcnt2_shifted
            extra_phase = UInt8(extra_phase_val & 0x3)
            @inbounds for b in 1:total_blocks
                d_mat[b, m] = d_mat[b, m] ⊻ s_pivot[b]
            end
            if update_phase
                d_phases[m] = (d_phases[m] + d_phases[current_col] + extra_phase) & UInt8(0x3)
            end
        end
    end
    return
end

function gpu_update_kernel!(d_mat::CuArray{T,2}, d_phases::CuVector{UInt8},
                           current_col::Integer, pivot_row::Integer, bit::T,
                           r::Integer, total_blocks::Integer, nblocks::Integer, 
                           update_phase::Bool) where {T}
    threads = 256
    blocks = max(cld(r, threads), MIN_BLOCKS)
    shmem_size = sizeof(T) * total_blocks
    @cuda threads=threads blocks=blocks shmem=shmem_size update_kernel!(
        d_mat, d_phases, current_col, pivot_row, bit, r,
        total_blocks, nblocks, update_phase)
    return
end

function gpu_canonicalize!(tableau::QuantumClifford.Tableau{<:CuArray{UInt8}, <:CuArray{T}}, phases::Bool) where {T<:Unsigned}
    d_mat_orig = tableau.xzs
    d_phases = tableau.phases
    nqubits = tableau.nqubits
    r = size(d_mat_orig, 2)    
    d_mat_flat = reinterpret(UInt64, reshape(d_mat_orig, (:,)))
    elements_per_row = size(d_mat_orig, 1) ÷ (sizeof(UInt64) ÷ sizeof(T))
    d_mat = reshape(d_mat_flat, (elements_per_row, size(d_mat_orig, 2)))    
    nblocks = cld(nqubits, 64)
    total_blk = 2 * nblocks
    current_col = 1    
    for j in 1:nqubits
        pivot_row = (j - 1) ÷ 64 + 1
        bit = UInt64(1) << ((j - 1) % 64)
        pivot = gpu_pivot_search!(d_mat, pivot_row, bit, current_col, r)
        if pivot == 0
            continue
        end
        if pivot != current_col
            if phases
                gpu_swap_columns_and_phases!(d_mat, d_phases, pivot, current_col)
            else
                gpu_swap_columns!(d_mat, pivot, current_col)
            end
        end
        gpu_update_kernel!(d_mat, d_phases, current_col, pivot_row, bit,
                          r, total_blk, nblocks, phases)
        current_col += 1
        if current_col > r
            break
        end
    end    
    for j in 1:nqubits
        pivot_row = nblocks + (j - 1) ÷ 64 + 1
        bit = UInt64(1) << ((j - 1) % 64)
        pivot = gpu_pivot_search!(d_mat, pivot_row, bit, current_col, r)
        if pivot == 0
            continue
        end
        if pivot != current_col
            if phases
                gpu_swap_columns_and_phases!(d_mat, d_phases, pivot, current_col)
            else
                gpu_swap_columns!(d_mat, pivot, current_col)
            end
        end
        gpu_update_kernel!(d_mat, d_phases, current_col, pivot_row, bit,
                          r, total_blk, nblocks, phases)
        current_col += 1
        if current_col > r
            break
        end
    end
    return tableau
end

function canonicalize!(stab::Stabilizer{<:QuantumClifford.Tableau{<:CuArray{UInt8}, <:CuArray{T}}}; phases::Bool=true) where {T<:Unsigned}
    gpu_canonicalize!(tab(stab), phases)
    CUDA.synchronize()
    return stab
end