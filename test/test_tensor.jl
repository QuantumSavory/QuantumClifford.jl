@testitem "Tensor product arity" begin
    using QuantumClifford
    using QuantumClifford: GeneralizedStabilizer
    using QuantumInterface: tensor

    @test !applicable(tensor)

    pauli = P"X"
    stabilizer = S"Z"
    destabilizer = MixedDestabilizer(stabilizer)
    clifford = CliffordOperator(sHadamard)
    register = Register(stabilizer, [true])
    generalized = GeneralizedStabilizer(stabilizer)

    @test tensor(pauli) == pauli
    @test tensor(stabilizer) == stabilizer
    @test tensor(destabilizer) == destabilizer
    @test tensor(clifford) == clifford
    @test tensor(register) == register
    @test tensor(pcT).paulis == collect(pcT.paulis)
    @test tensor(pcT).weights == collect(pcT.weights)

    @test tensor(register, stabilizer) == register ⊗ stabilizer
    @test tensor(generalized, stabilizer) == generalized ⊗ stabilizer
    @test tensor(pcT, pauli) == pcT ⊗ pauli

    @test tensor(pauli, P"Y", P"Z") == P"XYZ"
    @test tensor(stabilizer, S"X", S"Y") == stabilizer ⊗ S"X" ⊗ S"Y"
    @test tensor(clifford, tPhase, tHadamard) == clifford ⊗ tPhase ⊗ tHadamard
end
