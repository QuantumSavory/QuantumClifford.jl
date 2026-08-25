using TestItemRunner

@run_package_tests filter = ti -> (
    dirname(ti.filename) == (@__DIR__) &&
    (!(:bitpack in ti.tags) || (Sys.islinux() && Int === Int64))
)
