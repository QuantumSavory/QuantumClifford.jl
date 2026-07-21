Oscar_flag = false
Tesseract_flag = false
const JET_PROJECT = normpath(joinpath(@__DIR__, "projects", "jet"))
const test_args = isempty(ARGS) ? ["general"] : ARGS
const JET_flag = length(test_args) == 1 && startswith(only(test_args), "jet")

if JET_flag
    @info "Activating the dedicated JET test environment." project=JET_PROJECT
    using Pkg
    Pkg.activate(JET_PROJECT)
    Pkg.instantiate()
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
