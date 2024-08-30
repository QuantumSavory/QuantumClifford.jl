@testitem "throws" begin
    using QuantumClifford: rank, mul_left!, mul_right!
    using InteractiveUtils: subtypes

    @test_throws DimensionMismatch CliffordOperator(T"XXX ZZ_")

    @test_throws DimensionMismatch tCNOT*S"X"

    #@test_throws DomainError bigram(random_stabilizer(50), clip=false)

    @test_throws DomainError logdot(S"XX", S"XX ZZ")
    @test_throws DimensionMismatch logdot(S"X", S"XX ZZ")

    @test_throws BadDataStructure rank(S"X")
    @test_throws BadDataStructure rank(Destabilizer(S"X"))

    @test_throws DimensionMismatch mul_left!(P"X", P"XX")
    @test_throws DimensionMismatch mul_right!(P"X", P"XX")

    @test_throws ArgumentError StabMixture(S"XX")

    @test_throws ArgumentError UnitaryPauliChannel([P"X"], [1,2])
    @test_throws ArgumentError UnitaryPauliChannel([P"X",P"XX"], [1,2])

    @test_throws ArgumentError embed(10,2,P"XX")
    @test_throws ArgumentError embed(10,[2,3],P"X")

    struct A <: QuantumClifford.AbstractOperation end
    @test_throws ArgumentError applybranches(S"X",A())

    @test_throws BadDataStructure project!(Destabilizer(S"XX"), P"ZZ")

    @test_throws DimensionMismatch reset_qubits!(ghz(4), ghz(3), [1,2])
    @test_throws DimensionMismatch reset_qubits!(ghz(3), ghz(4), [1,2,3,4])
    @test_throws DimensionMismatch reset_qubits!(MixedStabilizer(ghz(4)), MixedStabilizer(ghz(3)), [1,2])
    @test_throws DimensionMismatch reset_qubits!(MixedStabilizer(ghz(3)), MixedStabilizer(ghz(4)), [1,2,3,4])
    @test_throws DimensionMismatch reset_qubits!(MixedDestabilizer(ghz(4)), MixedDestabilizer(ghz(3)), [1,2])
    @test_throws DimensionMismatch reset_qubits!(MixedDestabilizer(ghz(3)), MixedDestabilizer(ghz(4)), [1,2,3,4])

    #TODO broken in other ways @test_throws DomainError MixedDestabilizer(Destabilizer(S"XX"))

    @test_throws DomainError 2*P"X"

    @test_throws DimensionMismatch P"X" * S"XX"

    @test_throws ArgumentError one(typeof(T"X"), 1, basis=:U)

    for gt in subtypes(QuantumClifford.AbstractSingleQubitOperator)
        gt == SingleQubitOperator && continue
        @test_throws ArgumentError gt(0)
        @test_throws ArgumentError gt(-1)
    end

    for gt in subtypes(QuantumClifford.AbstractTwoQubitOperator)
        @test_throws ArgumentError gt(0,1)
        @test_throws ArgumentError gt(-1,1)
        @test_throws ArgumentError gt(2,2)
    end

    for m in [sMX,sMZ,sMY,sMRX,sMRZ,sMRY]
        @test_throws ArgumentError m(0)
        @test_throws ArgumentError m(-1)
        @test_throws ArgumentError m(0,1)
        @test_throws ArgumentError m(-1,0)
    end

    @test_throws DimensionMismatch Stabilizer(fill(0x0, 2), Matrix{Bool}([1 0 1;1 1 1; 1 0 1]), Matrix{Bool}([1 0 0;1 1 1;1 0 1]))
    @test_throws DimensionMismatch Stabilizer(fill(0x0, 2), Matrix{Bool}([1 0 1 1 0 0; 1 1 1 1 1 1; 1 0 1 1 0 1]))

end
