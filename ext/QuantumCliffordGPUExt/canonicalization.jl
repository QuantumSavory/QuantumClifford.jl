using CUDA
using LinearAlgebra
using QuantumClifford  
function swap_columns_kernel!(d_mat::CuDeviceArray{T,2}, 
                             d_phases::CuDeviceArray{P,1}, 
                             col1::Integer, col2::Integer, N::Integer,
                             update_phases::Bool) where {T, P}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx <= N
        temp = d_mat[idx, col1]
        d_mat[idx, col1] = d_mat[idx, col2]
        d_mat[idx, col2] = temp
    end
    if update_phases && idx == 1
        temp_phase = d_phases[col1]
        d_phases[col1] = d_phases[col2]
        d_phases[col2] = temp_phase
    end
    return
end

function gpu_swap_columns!(d_mat::CuArray{T,2}, d_phases::CuVector{P}, 
                          col1::Integer, col2::Integer, update_phases::Bool) where {T, P}
    N = size(d_mat, 1)
    threads = 256
    blocks = cld(N, threads)
    @cuda threads=threads blocks=blocks swap_columns_kernel!(
        d_mat, d_phases, col1, col2, N, update_phases)
    return
end

function pivot_search_kernel!(d_mat::CuDeviceArray{T,2}, pivot_row::Integer, bit::T, 
                             current_col::Integer, r::Integer, d_result) where {T}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x + current_col - 1
    
    if idx <= r
        if (d_mat[pivot_row, idx] & bit) != T(0)
            CUDA.atomic_min!(pointer(d_result), idx)
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
    d_result = CUDA.fill(r + 1, 1)
    threads = min(256, n_cols)
    blocks = cld(n_cols, threads)
    @cuda threads=threads blocks=blocks pivot_search_kernel!(
        d_mat, pivot_row, bit, current_col, r, d_result)
    result = Array(d_result)[1]
    return result > r ? 0 : Int(result)
end

function update_kernel!(d_mat::CuDeviceArray{T,2}, d_phases::CuDeviceArray{P,1},
                       current_col::Integer, pivot_row::Integer, bit::T,
                       r::Integer, total_blocks::Integer, nblocks::Integer, 
                       update_phase::Bool) where {T, P}
    s_pivot = @cuDynamicSharedMem(T, total_blocks)
    tid = threadIdx().x
    for idx in tid:blockDim().x:total_blocks
        s_pivot[idx] = d_mat[idx, current_col]
    end
    sync_threads()
    m = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if m <= r && m != current_col
        if (d_mat[pivot_row, m] & bit) != T(0)
            cnt1 = T(0)
            cnt2 = T(0)
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
            extra_phase = P(extra_phase_val & 0x3)
            @inbounds for b in 1:total_blocks
                d_mat[b, m] = d_mat[b, m] ⊻ s_pivot[b]
            end
            if update_phase
                d_phases[m] = (d_phases[m] + d_phases[current_col] + extra_phase) & P(0x3)
            end
        end
    end
    return
end

function gpu_update_kernel!(d_mat::CuArray{T,2}, d_phases::CuVector{P},
                           current_col::Integer, pivot_row::Integer, bit::T,
                           r::Integer, total_blocks::Integer, nblocks::Integer, 
                           update_phase::Bool) where {T, P}
    threads = 256
    blocks = cld(r, threads)
    shmem_size = sizeof(T) * total_blocks
    @cuda threads=threads blocks=blocks shmem=shmem_size update_kernel!(
        d_mat, d_phases, current_col, pivot_row, bit, r,
        total_blocks, nblocks, update_phase)
    return
end

function gpu_canonicalize!(tableau::QuantumClifford.Tableau{<:CuArray{P}, <:CuArray{T}}, phases::Bool) where {T, P}
    d_mat = tableau.xzs
    d_phases = tableau.phases
    nqubits = tableau.nqubits
    r = size(d_mat, 2)    
    nblocks = size(d_mat, 1) ÷ 2
    total_blk = 2 * nblocks
    bits_per_chunk = count_zeros(zero(T))
    current_col = 1
    for j in 1:nqubits
        pivot_row = (j - 1) ÷ bits_per_chunk + 1
        bit = T(1) << ((j - 1) % bits_per_chunk)
        pivot = gpu_pivot_search!(d_mat, pivot_row, bit, current_col, r)
        if pivot == 0
            continue
        end
        if pivot != current_col
            gpu_swap_columns!(d_mat, d_phases, pivot, current_col, phases)
        end
        gpu_update_kernel!(d_mat, d_phases, current_col, pivot_row, bit,
                          r, total_blk, nblocks, phases)
        current_col += 1
        if current_col > r
            break
        end
    end    
    for j in 1:nqubits
        pivot_row = nblocks + (j - 1) ÷ bits_per_chunk + 1
        bit = T(1) << ((j - 1) % bits_per_chunk)
        pivot = gpu_pivot_search!(d_mat, pivot_row, bit, current_col, r)
        if pivot == 0
            continue
        end
        if pivot != current_col
            gpu_swap_columns!(d_mat, d_phases, pivot, current_col, phases)
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

function canonicalize!(stab::Stabilizer{<:QuantumClifford.Tableau{<:CuArray{P}, <:CuArray{T}}}; phases::Bool=true) where {T, P}
    gpu_canonicalize!(tab(stab), phases)
    CUDA.synchronize()
    return stab
end