CUDA_flag = false
OpenCL_flag = false
ROCm_flag = false
Oscar_flag = false

if Sys.iswindows() || Sys.ARCH != :x86_64
    @info "Skipping GPU tests -- only supported on x86_64 *NIX platforms."
    @info "Skipping OpenCL tests -- only supported on x86_64 *NIX platforms."
    @info "Skipping Oscar tests -- only supported on x86_64 *NIX platforms."
else
    CUDA_flag = get(ENV, "CUDA_TEST", "") == "true"
    OpenCL_flag = true
    ROCm_flag = get(ENV, "ROCm_TEST", "") == "true"
    Oscar_flag = VERSION >= v"1.11"

    CUDA_flag && @info "Running with CUDA tests."
    OpenCL_flag && @info "Running with OpenCL tests."
    ROCm_flag && @info "Running with ROCm tests."
    !Oscar_flag && @info "Skipping Oscar tests -- not tested on Julia < 1.11"
    if !(ROCm_flag || CUDA_flag)
        @info "Skipping GPU tests -- must be explicitly enabled."
        @info "Environment must set ROCm_TEST=true xor CUDA_TEST=true."
    end
end

using Pkg
CUDA_flag && Pkg.add("CUDA")
OpenCL_flag && (Pkg.add("pocl_jll"); Pkg.add("OpenCL"))
ROCm_flag && Pkg.add("AMDGPU")
any((CUDA_flag, OpenCL_flag, ROCm_flag)) && Pkg.add("GPUArrays")
Oscar_flag && Pkg.add("Oscar")
using TestItemRunner
using QuantumClifford

# filter for the test
testfilter = ti -> begin
    exclude = Symbol[]

    if get(ENV, "JET_TEST", "") != "true"
        push!(exclude, :jet)
    else
        return :jet in ti.tags
    end

    if get(ENV, "ECC_TEST", "") != "true"
        push!(exclude, :ecc)
    else
        return :ecc in ti.tags
    end


    if !CUDA_flag
        push!(exclude, :cuda)
    else
        return :cuda in ti.tags
    end

    if !OpenCL_flag
        push!(exclude, :opencl)
    else
        return :opencl in ti.tags
    end

    if !ROCm_flag
        push!(exclude, :rocm)
    else
        return :rocm in ti.tags
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
