CUDA_flag = false
ROCm_flag = false
OpenCL_flag = false
Oscar_flag = false

if Sys.iswindows() || Sys.ARCH != :x86_64
    @info "Skipping Oscar tests -- only supported x86_64 *NIX platforms."
else
    Oscar_flag = VERSION >= v"1.11"
    !Oscar_flag && @info "Skipping Oscar tests -- not tested on Julia < 1.11"
end

if Sys.iswindows()
    @info "Skipping GPU/OpenCL tests -- only executed on *NIX platforms."
else
    CUDA_flag = get(ENV, "CUDA_TEST", "") == "true"
    ROCm_flag = get(ENV, "ROCm_TEST", "") == "true"
    OpenCL_flag = get(ENV, "OpenCL_TEST", "") == "true"

    CUDA_flag && @info "Running with CUDA tests."
    ROCm_flag && @info "Running with ROCm tests."
    OpenCL_flag && @info "Running with OpenCL tests."
    if !any((CUDA_flag, ROCm_flag, OpenCL_flag))
        @info "Skipping GPU/OpenCL tests -- must be explicitly enabled."
        @info "Environment must uniquely set [CUDA, ROCm, OpenCL]_TEST=true."
    end
end

using Pkg
CUDA_flag && Pkg.add("CUDA")
ROCm_flag && Pkg.add("AMDGPU")
OpenCL_flag && Pkg.add(["pocl_jll", "OpenCL"])
if any((CUDA_flag, ROCm_flag, OpenCL_flag))
    Pkg.add(
        ["Adapt", "Atomix", "GPUArraysCore", "GPUArrays", "KernelAbstractions"]
         )
end
Oscar_flag && Pkg.add("Oscar")
using TestItemRunner
using QuantumClifford

# filter for the test
testfilter = ti -> begin
    exclude = Symbol[]

    if get(ENV, "JET_TEST", "") == "true"
        return :jet in ti.tags
    else
        push!(exclude, :jet)
    end

    if get(ENV, "ECC_TEST_1", "") == "true"
        return (:ecc in ti.tags) && (:ecc_universal_checks in ti.tags)

    elseif get(ENV, "ECC_TEST_2", "") == "true"
        return (:ecc in ti.tags) && (:ecc_bespoke_checks in ti.tags)
    else
        push!(exclude, :ecc)
    end

    if CUDA_flag
        return :cuda in ti.tags
    else
        push!(exclude, :cuda)
    end

    if ROCm_flag
        return :rocm in ti.tags
    else
        push!(exclude, :rocm)
    end

    if OpenCL_flag
        return :opencl in ti.tags
    else
        push!(exclude, :opencl)
    end

    if !(VERSION >= v"1.10")
        push!(exclude, :doctests)
        push!(exclude, :aqua)
    end

    if !(Base.Sys.islinux() & (Int===Int64))
        push!(exclude, :bitpack)
    end

    return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
