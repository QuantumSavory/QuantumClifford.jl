@testitem "JET checks" tags=[:jet] begin
    using JET
    using Test
    using QECCore

    JET.test_package(QECCore, target_modules=[QECCore])
end
