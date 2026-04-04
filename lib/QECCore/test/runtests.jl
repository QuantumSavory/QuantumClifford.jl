using QECCore
using TestItemRunner

Oscar_flag = false
Tesseract_flag = false
JET_flag = false

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

if get(ENV, "JET_TEST", "") == "true"
    JET_flag = true
else
    @info "Skipping JET tests -- must be explicitly enabled."
end

using Pkg
Oscar_flag && Pkg.add("Oscar")
Tesseract_flag && Pkg.add("PyTesseractDecoder")
JET_flag && Pkg.add("JET")

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
