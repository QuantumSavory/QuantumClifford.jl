using CUDA
using BenchmarkTools
using LinearAlgebra
using QuantumClifford  

const MIN_BLOCKS = 36

###########################################################################
# Custom Linear Copy Kernel
###########################################################################
# This kernel copies elements from a source device array (src) to a destination device array (dest).
function linear_copy_kernel!(dest, src, N)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    while i <= N
        dest[i] = src[i]
        i += stride
    end
    return
end

# A helper function that first converts the host array to a device array and then uses our custom kernel.
function custom_gpu_linear_copy(arr::Array{T}) where T
    N = length(arr)
    # Convert host array to device array (this call is usually optimized with cudaMemcpy)
    d_arr_src = CuArray(arr)
    d_arr_dest = CUDA.zeros(T, N)
    threads = 256  # a multiple of 32
    blocks = MIN_BLOCKS  
    @cuda threads=threads blocks=blocks linear_copy_kernel!(d_arr_dest, d_arr_src, N)
    return d_arr_dest
end

###########################################################################
# 1. Data Structure: GpuTableau
###########################################################################
struct GpuTableau
    d_mat::CuArray{UInt64,2}    # Bit-packed matrix on GPU.
                                # Rows 1:nblocks hold the X–part;
                                # rows nblocks+1:2*nblocks hold the Z–part.
    d_phases::CuVector{UInt8}   # Phases stored on GPU.
    nqubits::Int                # Number of qubits.
    nblocks::Int                # cld(nqubits, 64)
    r::Int                      # Number of generators (columns)
end

###########################################################################
# 2. GPU Kernels and Device Functions
###########################################################################

# --- Merged Swap Kernel (columns and phase vector) ---
function swap_columns_and_phases_kernel!(d_mat, d_phases, col1::Int32, col2::Int32, N::Int32)
    idx = (Int32(blockIdx().x) - 1) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx <= N
        temp = d_mat[idx, col1]
        d_mat[idx, col1] = d_mat[idx, col2]
        d_mat[idx, col2] = temp
    end
    # Only one thread swaps the phases.
    if idx == 1
        temp_phase = d_phases[col1]
        d_phases[col1] = d_phases[col2]
        d_phases[col2] = temp_phase
    end
    return
end

function gpu_swap_columns_and_phases!(d_mat, d_phases, col1::Int, col2::Int)
    N = size(d_mat, 1)
    threads = 256  # must be multiple of 32, tests revealed 256 is optimal
    blocks = max(cld(N, threads), MIN_BLOCKS)
    @cuda threads=threads blocks=blocks swap_columns_and_phases_kernel!(d_mat, d_phases,
        Int32(col1), Int32(col2), Int32(N))
    return
end

# --- Original Swap Kernel (matrix-only) ---
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

# --- Popcount and Compute Phase (CPU–side versions for reference) ---
@inline function gpu_popc_32(x::UInt32)
    count = UInt32(0)
    @inbounds for i in 0:31
        count += (x >> i) & UInt32(1)
    end
    return count
end

@inline function gpu_popc(x::UInt64)
    lo = gpu_popc_32(UInt32(x))
    hi = gpu_popc_32(UInt32(x >> 32))
    return lo + hi
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
    rcnt1 = gpu_popc(cnt1)
    rcnt2 = gpu_popc(cnt2)
    return UInt8((rcnt1 ⊻ (rcnt2 << 1)) & UInt32(0x3))
end

# --- Pivot Search Kernel ---
function pivot_search_kernel!(d_mat, pivot_row::Int32, bit::UInt64, current_col::Int32, r::Int32, d_result)
    result = Int32(0)
    m = current_col
    while m <= r
        if (d_mat[pivot_row, m] & bit) != UInt64(0)
            result = m
            break
        end
        m += Int32(1)
    end
    d_result[1] = result
    return
end

function gpu_pivot_search!(d_mat, pivot_row::Int, bit::UInt64, current_col::Int, r::Int)::Int
    row_segment = Array(view(d_mat, pivot_row, current_col:r))
    pivot_offset = findfirst(x -> (x & bit) != UInt64(0), row_segment)
    if pivot_offset === nothing
        return 0
    else
        return pivot_offset + current_col - 1
    end
end

# --- Update Kernel with Shared Memory Caching ---
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
            lo1 = UInt32(cnt1 & 0xffffffff)
            hi1 = UInt32((cnt1 >> 32) & 0xffffffff)
            count_lo1 = UInt32(0)
            count_hi1 = UInt32(0)
            @inbounds for i in 0:31
                count_lo1 += (lo1 >> i) & UInt32(1)
                count_hi1 += (hi1 >> i) & UInt32(1)
            end
            rcnt1 = count_lo1 + count_hi1
            lo2 = UInt32(cnt2 & 0xffffffff)
            hi2 = UInt32((cnt2 >> 32) & 0xffffffff)
            count_lo2 = UInt32(0)
            count_hi2 = UInt32(0)
            @inbounds for i in 0:31
                count_lo2 += (lo2 >> i) & UInt32(1)
                count_hi2 += (hi2 >> i) & UInt32(1)
            end
            rcnt2 = count_lo2 + count_hi2
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

###########################################################################
# 3. GPU Canonicalization Routine
###########################################################################
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

###########################################################################
# 4. GPU-Resident Canonicalization Wrapper
###########################################################################
function canonicalize_gpu!(stab::QuantumClifford.Stabilizer; phases::Bool=true)
    d_mat    = stab.tab.xzs
    d_phases = stab.tab.phases
    nqubits  = stab.tab.nqubits
    r        = size(d_mat, 2)
    nblocks  = cld(nqubits, 64)
    gt = GpuTableau(d_mat, d_phases, nqubits, nblocks, r)
    gpu_canonicalize!(gt, phases)
    CUDA.synchronize()
    return stab
end

###########################################################################
# 5. Utilities: Converting between CPU and GPU Stabilizers
###########################################################################
function cpu_to_gpu_stabilizer(stab::QuantumClifford.Stabilizer)
    nqubits = stab.tab.nqubits
    phases  = stab.tab.phases
    B       = stab.tab.xzs  # may be Boolean or already bit-packed?
    bitpacked = (eltype(B) == Bool) ? bool_to_bitpacked(B, nqubits) : B
    d_mat = reshape(custom_gpu_linear_copy(vec(bitpacked)), size(bitpacked))
    d_phases = custom_gpu_linear_copy(phases)
    gpu_tab  = QuantumClifford.Tableau(d_phases, nqubits, d_mat)
    return QuantumClifford.Stabilizer(gpu_tab)
end

function get_stabilizer_matrix(gt::GpuTableau)
    r = gt.r
    nqubits = gt.nqubits
    nblocks = gt.nblocks
    M = Array(gt.d_mat)
    S = Array{Bool}(undef, r, 2 * nqubits)
    for i in 1:r
        for j in 1:nqubits
            block = div(j - 1, 64) + 1
            bit_pos = (j - 1) % 64
            S[i, j]         = ((M[block, i] >> bit_pos) & 1) == 1
            S[i, nqubits+j] = ((M[nblocks+block, i] >> bit_pos) & 1) == 1
        end
    end
    return S
end

function gpu_to_cpu_stabilizer(gpu_stab::QuantumClifford.Stabilizer)
    t = gpu_stab.tab
    nqubits = t.nqubits
    nblocks = cld(nqubits, 64)
    gt = GpuTableau(t.xzs, t.phases, nqubits, nblocks, size(t.xzs, 2))
    bool_mat = get_stabilizer_matrix(gt)
    xs = bool_mat[:, 1:nqubits]
    zs = bool_mat[:, nqubits+1:end]
    cpu_phases = Array(t.phases)
    new_tab = QuantumClifford.Tableau(cpu_phases, xs, zs)
    return QuantumClifford.Stabilizer(new_tab)
end

function gpu_stabilizer_copy(stab::QuantumClifford.Stabilizer)
    new_xzs = CUDA.copy(stab.tab.xzs)
    new_phases = CUDA.copy(stab.tab.phases)
    new_tab = QuantumClifford.Tableau(new_phases, stab.tab.nqubits, new_xzs)
    return QuantumClifford.Stabilizer(new_tab)
end


###########################################################################
# 7. Example Run and Verification
###########################################################################
if abspath(PROGRAM_FILE) == @__FILE__
    #String example
cpu_stab_example = S"""
+ Z_XY__YZXXYYZZ_Y__XZ_XZZ_X_Z__
+ _YZZXX_YXYZXYZ_ZZXZYXZYZ_ZZYZZ
- _YYZYX_YYXXZ__X__Y_ZYYYXX_XYYX
- _XZXYYYYY_Y_YY_ZY_XXZ_ZYXY_ZZY
- _YX_XZXXX_YY_ZXXYYYXZ_Y_YZ_YZY
+ Y__YZZXZXZ_Y_YXXXZ_Y_YZZZYZ_X_
+ YYYZZYXYYZYXXZZYXY__XXZZXYZ_ZX
+ YYZXYZ_XY___ZYZZZXXZZZZYYY_ZXY
- XYYX_Y_YY_ZX__Y_ZX__YYZXZ__Y__
- _ZYYXYYY_XZZY_ZZY_XY_X_Z_X_XY_
- XZ__XZYZYZYYX__X_ZZYXXXZZXY__Z
- _ZXZXZY___ZYYZYXXX_Y_XZ__ZXZZZ
+ ZYZZXX_ZXYZYXYZXYZYX_ZXZZZYZ__
+ X_YXXY_XXYZXYYZXYZYX_ZZY__ZZ_X
+ YZYYXXZYXZZ_Y_XXX_YYZ__XYYZZ_Y
- Z___YZZZXXXYYY_ZYXY_YZXY__XYY_
- ZXYYX_ZZY_Z_XZX_XXZYZZZZY____X
+ __XZZYXYZX_X_YZ_XZZZYXY_ZZZYY_
- YZXXYZ__Z_YZZXY_ZYZZZ_YZYYZ_Z_
+ Y_Y_XZY__ZZ_XY_Y_X__YZ_ZY__XXX
- X_Z_YX_ZXZYYXX_YY_ZXY__XZZZXZY
+ ZZXZXYZYXZYZZZYYZYZYX__X_YYX_X
- YZ__XZYYZXZXZ__XY_XXX_Y_YYZYXY
- Y__Z_XZ_Y_ZZXXZYZXYZYY_ZZZ_YX_
- ZY_YYXZXXXZZYXZ__X_ZYYXYZZZZ_Z
- ZZX_Y_XZXYZXZZ_XZ__YZXZX_YYZYZ
- __ZZZ__YYZXZZXXYYZYYXYZ_YXYX__
+ ZXYYZ_Z_XZ__YZ_YZYZYX__ZYZ_YX_
- ZYZZYYXYZYZ_ZZ__ZX__ZZ_ZXZ_YZX
+ Z_YYY_ZZ_X_Z_Y__XX_YYXXX__XXXX
"""
gpu_stab_example = cpu_to_gpu_stabilizer(cpu_stab_example)
canonicalize_gpu!(gpu_stab_example; phases=true)
println(typeof(gpu_stab_example))
cpu_stab_converted = gpu_to_cpu_stabilizer(gpu_stab_example)
println("Original CPU Stabilizer:")
println(cpu_stab_example)
println("\nCanonicalized Stabilizer (converted back from GPU):")
println(cpu_stab_converted)
println(typeof(cpu_stab_converted))

#using random_stabilizer
cpu_stabil = random_stabilizer(4)
gpu_stabil = cpu_to_gpu_stabilizer(cpu_stabil)
canonicalize_gpu!(gpu_stabil; phases = true)
cpu_stab_con = gpu_to_cpu_stabilizer(gpu_stabil)
println("Original CPU Stabilizer:")
println(cpu_stabil)
println("\nCanonicalized Stabilizer (converted back from GPU):")
println(cpu_stab_con)

end

