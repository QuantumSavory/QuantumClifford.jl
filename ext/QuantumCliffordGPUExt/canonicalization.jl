using CUDA
using LinearAlgebra
using QuantumClifford  

const MIN_BLOCKS = 36
const GPUTableau = QuantumClifford.Tableau{CuArray{UInt8, 1, CUDA.Mem.DeviceBuffer}, CuArray{UInt32, 2, CUDA.Mem.DeviceBuffer}}
const GPUStabilizer = QuantumClifford.Stabilizer{GPUTableau}

struct GpuTableau
    d_mat::CuArray{UInt64,2}         
    d_phases::CuVector{UInt8}   
    nqubits::Int                
    nblocks::Int                
    r::Int                      
end

function swap_columns_and_phases_kernel!(d_mat, d_phases, col1::Int32, col2::Int32, N::Int32)
    idx = (Int32(blockIdx().x) - 1) * Int32(blockDim().x) + Int32(threadIdx().x)
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

function gpu_swap_columns_and_phases!(d_mat, d_phases, col1::Int, col2::Int)
    N = size(d_mat, 1)
    threads = 256 
    blocks = max(cld(N, threads), MIN_BLOCKS)
    @cuda threads=threads blocks=blocks swap_columns_and_phases_kernel!(d_mat, d_phases,
        Int32(col1), Int32(col2), Int32(N))
    return
end

function swap_columns_kernel!(d_mat, col1::Int32, col2::Int32, N::Int32)
    idx = (Int32(blockIdx().x) - 1) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx <= N
        temp = d_mat[idx, col1]
        d_mat[idx, col1] = d_mat[idx, col2]
        d_mat[idx, col2] = temp
    end
    return
end

function gpu_swap_columns!(d_mat, col1::Int, col2::Int)
    N = size(d_mat, 1)
    threads = 256
    blocks = max(cld(N, threads), MIN_BLOCKS)
    @cuda threads=threads blocks=blocks swap_columns_kernel!(d_mat, Int32(col1), Int32(col2), Int32(N))
    return
end


@inline function gpu_compute_phase(d_mat, pivot_col::Int32, target_col::Int32, nblocks::Int32)::UInt8
    cnt1 = UInt64(0)
    cnt2 = UInt64(0)
    @inbounds for i in Int32(1):nblocks
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
    return UInt8((rcnt1 ⊻ (rcnt2 << 1)) & UInt32(0x3))
end

function pivot_search_kernel!(d_mat, pivot_row::Int32, bit::UInt64, current_col::Int32, r::Int32, d_result)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x + current_col - 1
    
    if idx <= r
        if (d_mat[pivot_row, idx] & bit) != UInt64(0)
            CUDA.atomic_min!(pointer(d_result), Int32(idx))
        end
    end
    return
end

function gpu_pivot_search!(d_mat, pivot_row::Int, bit::UInt64, current_col::Int, r::Int)::Int
    n_cols = r - current_col + 1
    if n_cols <= 0
        return 0
    end
    
    d_result = CUDA.fill(Int32(r + 1), 1)
    
    threads = min(256, n_cols)
    blocks = cld(n_cols, threads)
    
    @cuda threads=threads blocks=blocks pivot_search_kernel!(
        d_mat, Int32(pivot_row), bit, Int32(current_col), Int32(r), d_result
    )
    
    result = Array(d_result)[1]
    return result > r ? 0 : Int(result)
end

function update_kernel!(d_mat, d_phases, current_col::Int32, pivot_row::Int32, bit::UInt64,
                        r::Int32, total_blocks::Int32, nblocks::Int32, update_phase::Bool)
    s_pivot = @cuDynamicSharedMem(UInt64, total_blocks)
    tid = threadIdx().x
    for idx in tid: blockDim().x: total_blocks
        s_pivot[idx] = d_mat[idx, current_col]
    end
    sync_threads()
    m = (Int32(blockIdx().x) - 1) * blockDim().x + Int32(threadIdx().x)
    if m <= r && m != current_col
        if (d_mat[pivot_row, m] & bit) != UInt64(0)
            cnt1 = UInt64(0)
            cnt2 = UInt64(0)
            @inbounds for i in Int32(1):nblocks
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
            mask = UInt32(0x3)
            extra_phase = UInt8(extra_phase_val & mask)
            @inbounds for b in Int32(1):total_blocks
                d_mat[b, m] = d_mat[b, m] ⊻ s_pivot[b]
            end
            if update_phase
                d_phases[m] = (d_phases[m] + d_phases[current_col] + extra_phase) & UInt8(0x3)
            end
        end
    end
    return
end

function gpu_update_kernel!(d_mat, d_phases, current_col::Int, pivot_row::Int, bit::UInt64,
                              r::Int, total_blocks::Int, nblocks::Int, update_phase::Bool)
    threads = 256
    blocks = max(cld(r, threads), MIN_BLOCKS)
    shmem_size = sizeof(UInt64) * total_blocks
    @cuda threads=threads blocks=blocks shmem=shmem_size update_kernel!(d_mat, d_phases,
        Int32(current_col), Int32(pivot_row), bit, Int32(r),
        Int32(total_blocks), Int32(nblocks), update_phase)
    return
end

function gpu_canonicalize!(gt::GpuTableau, phases::Bool)
    d_mat     = gt.d_mat
    d_phases  = gt.d_phases
    r_host    = gt.r
    nqubits_host = gt.nqubits
    nblocks_host = gt.nblocks
    r        = Int32(r_host)
    nqubits  = Int32(nqubits_host)
    nblocks  = Int32(nblocks_host)
    total_blk = Int32(2) * nblocks
    current_col = Int32(1)
    for j in Int32(1):nqubits
        pivot_row = (j - Int32(1)) ÷ Int32(64) + Int32(1)
        bit = UInt64(1) << ((j - Int32(1)) % Int32(64))
        pivot = gpu_pivot_search!(d_mat, Int(pivot_row), bit, Int(current_col), Int(r))
        if pivot == 0
            continue
        end
        if pivot != current_col
            if phases
                gpu_swap_columns_and_phases!(d_mat, d_phases, Int(pivot), Int(current_col))
            else
                gpu_swap_columns!(d_mat, Int(pivot), Int(current_col))
            end
        end
        gpu_update_kernel!(d_mat, d_phases, Int(current_col), Int(pivot_row), bit,
                           Int(r), Int(total_blk), Int(nblocks), phases)
        current_col += Int32(1)
        if current_col > r
            break
        end
    end
    for j in Int32(1):nqubits
        pivot_row = nblocks + (j - Int32(1)) ÷ Int32(64) + Int32(1)
        bit = UInt64(1) << ((j - Int32(1)) % Int32(64))
        pivot = gpu_pivot_search!(d_mat, Int(pivot_row), bit, Int(current_col), Int(r))
        if pivot == 0
            continue
        end
        if pivot != current_col
            if phases
                gpu_swap_columns_and_phases!(d_mat, d_phases, Int(pivot), Int(current_col))
            else
                gpu_swap_columns!(d_mat, Int(pivot), Int(current_col))
            end
        end
        gpu_update_kernel!(d_mat, d_phases, Int(current_col), Int(pivot_row), bit,
                           Int(r), Int(total_blk), Int(nblocks), phases)
        current_col += Int32(1)
        if current_col > r
            break
        end
    end
    return gt
end

function canonicalize!(stab::GPUStabilizer; phases::Bool=true)
    d_mat_32 = stab.tab.xzs
    d_phases = stab.tab.phases
    nqubits  = stab.tab.nqubits
    r        = size(d_mat_32, 2)
    nblocks  = cld(nqubits, 64)
    
    d_mat_64 = reinterpret(UInt64, reshape(d_mat_32, (:,)))
    d_mat_64 = reshape(d_mat_64, (size(d_mat_32, 1) ÷ 2, size(d_mat_32, 2)))
    
    gt = GpuTableau(d_mat_64, d_phases, nqubits, nblocks, r)
    gpu_canonicalize!(gt, phases)
    CUDA.synchronize()
    return stab
end