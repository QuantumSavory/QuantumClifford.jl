using TestItemRunner

include("TestItemConfiguration.jl")
using .TestItemConfiguration: platform_supported

# Bit-packing runs only on 64-bit Linux; all gates are centralized here.
@run_package_tests filter = ti -> (
    dirname(ti.filename) == (@__DIR__) &&
    platform_supported(ti.tags)
)
