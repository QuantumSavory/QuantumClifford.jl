@testitem "Apply Right" begin
    using InteractiveUtils
    using QuantumClifford: AbstractCliffordOperator, CliffordOperator
    
    test_sizes = [1,2,63,64,65,127,128,129,511,512,513]

    function CliffordOperator(pauli::PauliOperator)
        n = nqubits(pauli)
        res = one(CliffordOperator, n)

        for i in 1:2n
            if comm(pauli, tab(res), i) == 0x1
                phases(res)[i] ⊻= 0x02
            end
        end
        return res
    end
    
    @testset "CliffordOperator constructor from PauliOperator" begin
        for n in test_sizes
            l = random_clifford(n)
            pauli = random_pauli(n)
            @test isequal(apply!(copy(l), pauli; phases=true), apply!(l, CliffordOperator(pauli); phases=true))
        end
    end
    
    # SLOW version of apply_right! for testing
    function apply_right_slow!(l::CliffordOperator, r::CliffordOperator; phases=true)
        apply!(r, l; phases)
    end
    function apply_right_slow!(l::CliffordOperator, r::PauliOperator; phases=true)
        apply!(CliffordOperator(r), l; phases)
    end
    function apply_right_slow!(l::CliffordOperator, r::AbstractCliffordOperator; phases=true)
        apply!(CliffordOperator(r, nqubits(l)), l; phases)
    end

    @testset "Apply Right dense operators" begin
        for q in test_sizes
            l = random_clifford(q)
            cliff = random_clifford(q)
            pauli = random_pauli(q)

            @test isequal(apply_right!(copy(l), cliff), apply_right_slow!(l, cliff))
            @inferred apply_right!(copy(l), cliff)

            @test isequal(apply_right!(copy(l), pauli), apply_right_slow!(l, pauli))
            @inferred apply_right!(copy(l), pauli)
        end
    end

    @testset "Apply Right single-qubit" begin
        for q in test_sizes
            l = random_clifford(q)
            q1 = rand(1:q)

            for gate in subtypes(AbstractSingleQubitOperator)
                gate == SingleQubitOperator && continue

                @test isequal(apply_right!(copy(l), gate(q1)), apply_right_slow!(l, gate(q1)))
                @inferred apply_right!(copy(l), gate(q1))
            end
        end
    end

    @testset "Apply Right two-qubit" begin
        for q in test_sizes
            q < 2 && continue

            l = random_clifford(q)
            q1 = rand(1:q); q2 = rand(setdiff(1:q, [q1]))

            for gate in subtypes(AbstractTwoQubitOperator)
                @test isequal(apply_right!(copy(l), gate(q1, q2)), apply_right_slow!(l, gate(q1, q2)))
                @inferred apply_right!(copy(l), gate(q1, q2))
            end
        end
    end
end