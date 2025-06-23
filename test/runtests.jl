using Pkg
if Sys.iswindows() || Sys.ARCH != :x86_64
    @info "skipping Oscar tests (they currently do not run on Windows OS or ARM CPU)"
    @info "skipping GPU tests (set GPU_TESTS=true to test GPU (on non-Windows))"
elseif get(ENV, "GPU_TESTS", "") == "true"
    @info "running with GPU tests"
    Pkg.add("CUDA")
elseif VERSION < v"1.11"
    @info "skipping Oscar tests (not tested on Julia <1.11)"
    @info "skipping GPU tests (set GPU_TESTS=true to test GPU)"
else
    @info "skipping GPU tests (set GPU_TESTS=true to test GPU)"
    Pkg.add("Oscar")
end

using TestItemRunner
using QuantumClifford

# filter for the test
testfilter = ti -> begin
    exclude = Symbol[]

    if get(ENV, "JET_TEST", "")!="true"
        push!(exclude, :jet)
    else
        return :jet in ti.tags
    end

    if !(VERSION >= v"1.10")
        push!(exclude, :doctests)
        push!(exclude, :aqua)
    end

    if get(ENV, "GPU_TESTS", "")!="true"
        push!(exclude, :gpu)
    end

    if !(Base.Sys.islinux() & (Int===Int64))
        push!(exclude, :bitpack)
    end

    return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
