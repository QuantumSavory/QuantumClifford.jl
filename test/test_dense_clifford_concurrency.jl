@testitem "Dense Clifford application supports concurrent reuse" begin
    using Random: Xoshiro
    using QuantumClifford

    gate = random_clifford(Xoshiro(11), 5)
    inputs = [random_stabilizer(Xoshiro(seed), 5) for seed in 21:36]
    expected = [apply!(copy(state), copy(gate)) for state in inputs]
    results = similar(expected)

    @sync for i in eachindex(inputs)
        Threads.@spawn results[i] = apply!(copy(inputs[i]), gate)
    end
    @test results == expected

    indexed_gate = random_clifford(Xoshiro(41), 2)
    indices = [2, 5]
    expected = [apply!(copy(state), indices, copy(indexed_gate)) for state in inputs]
    results = similar(expected)
    @sync for i in eachindex(inputs)
        Threads.@spawn results[i] = apply!(copy(inputs[i]), indices, indexed_gate)
    end
    @test results == expected
end
