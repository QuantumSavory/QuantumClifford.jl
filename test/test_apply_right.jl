@testitem "Apply Right" begin
    using InteractiveUtils
    using QuantumClifford: AbstractCliffordOperator
    
    # SLOW version of apply_right! for testing
    function apply_right_slow!(l::CliffordOperator, r::AbstractCliffordOperator; phases=true)
        apply!(CliffordOperator(r, nqubits(l)), l; phases=phases)
    end
    
    q = 64
    shots = 16

    @testset "Apply Right single-qubit" begin
        for _ in 1:shots
            l = random_clifford(q)
            q1 = rand(1:q)

            for gate in subtypes(AbstractSingleQubitOperator)
                gate == SingleQubitOperator && continue
                @test isequal(apply_right!(copy(l), gate(q1)), apply_right_slow!(l, gate(q1)))
            end
        end
    end

    @testset "Apply Right two-qubit" begin
        for _ in 1:shots
            l = random_clifford(q)
            q1 = rand(1:q); q2 = rand(setdiff(1:q, [q1]))

            for gate in subtypes(AbstractTwoQubitOperator)
                @test isequal(apply_right!(copy(l), gate(q1, q2)), apply_right_slow!(l, gate(q1, q2)))
            end
        end
    end
end