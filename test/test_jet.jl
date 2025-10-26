@testitem "JET checks" tags=[:jet] begin
    using JET
    using Test
    using QuantumClifford
    using QuantumInterface

    JET.test_package(QuantumClifford, target_modules=[QuantumClifford, QuantumInterface])
end
