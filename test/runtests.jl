AMDGPU_flag = false
CUDA_flag = false
Oscar_flag = false

if Sys.iswindows() || Sys.ARCH != :x86_64
    @info "Skipping GPU tests -- only supported on x86_64 *NIX platforms."
    @info "Skipping Oscar tests -- only supported on x86_64 *NIX platforms."
else
    AMDGPU_flag = get(ENV, "AMDGPU_TEST", "") == "true"
    CUDA_flag = get(ENV, "CUDA_TEST", "") == "true"
    Oscar_flag = VERSION >= v"1.11"

    AMDGPU_flag && @info "Running with AMDGPU tests."
    CUDA_flag && @info "Running with CUDA tests."
    !Oscar_flag && @info "Skipping Oscar tests -- not tested on Julia < 1.11"
    if !(AMDGPU_flag || CUDA_flag)
        @info "Skipping GPU tests -- must be explicitly enabled."
        @info "Environment must set AMDGPU_TEST xor CUDA_TEST to \"true\"."
    end
end

using Pkg
AMDGPU_flag && Pkg.add("AMDGPU")
CUDA_flag && Pkg.add("CUDA")
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

    if !AMDGPU_flag
        push!(exclude, :amdgpu)
    else
        return :amdgpu in ti.tags
    end

    if !CUDA_flag
        push!(exclude, :cuda)
    else
        return :cuda in ti.tags
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
