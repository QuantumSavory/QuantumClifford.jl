using Pkg
if Sys.iswindows()
    @info "skipping Oscar tests (they currently do not run on Windows)"
    @info "skipping GPU tests (set GPU_TESTS=true to test GPU (on non-Windows))"
elseif get(ENV, "GPU_TESTS", "") == "true"
    Pkg.activate("gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @info "running with GPU tests"
elseif VERSION < v"1.11"
    @info "skipping Oscar tests (not tested on Julia <1.11)"
    @info "skipping GPU tests (set GPU_TESTS=true to test GPU)"
else
    @info "skipping GPU tests (set GPU_TESTS=true to test GPU)"
    Pkg.activate("default")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

using TestItemRunner
using QuantumClifford

# filter for the test
testfilter = ti -> begin
    exclude = Symbol[]
    if get(ENV, "JET_TEST", "")!="true"
        push!(exclude, :jet)
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
