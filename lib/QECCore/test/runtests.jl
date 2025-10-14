using QECCore
using TestItemRunner

Oscar_flag = false

if Sys.iswindows() || Sys.ARCH != :x86_64
    @info "Skipping Oscar tests -- only supported x86_64 *NIX platforms."
else
    Oscar_flag = VERSION >= v"1.11"
    !Oscar_flag && @info "Skipping Oscar tests -- not tested on Julia < 1.11"
end

using Pkg
Oscar_flag && Pkg.add("Oscar")

@run_package_tests
