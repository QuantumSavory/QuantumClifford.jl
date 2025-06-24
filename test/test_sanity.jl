@testitem "Simple sanity checks" begin
    using Random
    using InteractiveUtils
    @testset "Apply small symbolics" begin
        for gate_type in subtypes(QuantumClifford.AbstractTwoQubitOperator)
            @test apply!(S"II II", gate_type(1,2)) == S"II II"
            s1 = random_stabilizer(4)
            s2 = s1[randperm(4)]
            apply!(s1, gate_type(3,2))
            apply!(s2, gate_type(3,2))
            @test canonicalize!(s1) == canonicalize!(s2)
        end
        for gate_type in subtypes(QuantumClifford.AbstractSingleQubitOperator)
            gate_type === SingleQubitOperator && continue
            @test apply!(S"I", gate_type(1)) == S"I"
            s1 = random_stabilizer(4)
            s2 = s1[randperm(4)]
            apply!(s1, gate_type(2))
            apply!(s2, gate_type(2))
            @test canonicalize!(s1) == canonicalize!(s2)
            gate(x) = SingleQubitOperator(gate_type(x))
            @test apply!(S"I", gate(1)) == S"I"
            s1 = random_stabilizer(4)
            s2 = s1[randperm(4)]
            apply!(s1, gate(2))
            apply!(s2, gate(2))
            @test canonicalize!(s1) == canonicalize!(s2)
        end
    end

    @testset "Apply small Cliffords" begin
        for gate_type in (random_pauli, random_clifford)
            @test apply!(S"II II", gate_type(2)) == S"II II"
            s1 = random_stabilizer(4)
            s2 = s1[randperm(4)]
            gate = gate_type(4)
            apply!(s1, gate)
            apply!(s2, gate)
            @test canonicalize!(s1) == canonicalize!(s2)
        end
    end
end
