@testitem "All-to-all circuits use the supplied RNG" begin
    using Random
    using QuantumClifford

    signature(circuit) = [(copy(g.indices), copy(tab(g.cliff))) for g in circuit]

    Random.seed!(Random.GLOBAL_RNG, 71)
    first_circuit = random_all_to_all_clifford_circuit(Xoshiro(72), 8, 20)
    Random.seed!(Random.GLOBAL_RNG, 73)
    second_circuit = random_all_to_all_clifford_circuit(Xoshiro(72), 8, 20)

    @test signature(first_circuit) == signature(second_circuit)
end
