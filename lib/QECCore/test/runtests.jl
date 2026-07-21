Oscar_flag = false
Tesseract_flag = false
JET_flag = ARGS == ["jet"] || get(ENV, "JET_TEST", "") == "true"

if JET_flag
    @info "Running JET tests in their dedicated test environment."
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "projects", "jet"))
    Pkg.instantiate()
else
    @info "Skipping JET tests -- pass `test_args=[\"jet\"]` to Pkg.test to enable them."
end

using QECCore
using TestItemRunner

if Sys.iswindows() || Sys.ARCH != :x86_64
    @info "Skipping Oscar tests -- only supported x86_64 *NIX platforms."
else
    Oscar_flag = VERSION >= v"1.11"
    !Oscar_flag && @info "Skipping Oscar tests -- not tested on Julia < 1.11"
    Tesseract_flag = true
end

if Sys.iswindows()
    @info "Skipping Tesseract tests -- only supported *NIX platforms."
else
    Tesseract_flag = true
end

using Pkg
!JET_flag && Oscar_flag && Pkg.add("Oscar")
!JET_flag && Tesseract_flag && Pkg.add("PyTesseractDecoder")

# filter for the test
testfilter = ti -> begin
    exclude = Symbol[]

    if JET_flag
        return :jet in ti.tags
    else
        push!(exclude, :jet)
    end

    if !(VERSION >= v"1.10")
        push!(exclude, :doctests)
        push!(exclude, :aqua)
    end

    return all(!in(exclude), ti.tags)
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@run_package_tests filter=testfilter
