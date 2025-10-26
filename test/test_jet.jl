@testitem "JET checks" tags=[:jet] begin
    using JET
    using Test
    using QuantumClifford

    JET.test_package(QuantumClifford)
end
