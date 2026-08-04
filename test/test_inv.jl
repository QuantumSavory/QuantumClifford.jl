@testitem "Apply Inv" begin
    using InteractiveUtils

    @testset "Apply Inv" begin
        stabilizers = [S" I", S" X", S" Y", S" Z",
                       S"-I", S"-X", S"-Y", S"-Z",]
        for gate in subtypes(AbstractSingleQubitOperator)
            gate == SingleQubitOperator && continue
            for stab in stabilizers
                @test apply_inv!(stab, gate(1)) == apply!(stab, inv(CliffordOperator(gate(1), 1)))
            end
        end
    end

    @testset "Apply Inv two-qubit" begin
        stabilizers = [S" II", S" IX", S" IY", S" IZ",
                       S"-II", S"-IX", S"-IY", S"-IZ",
                       S" XI", S" XX", S" XY", S" XZ",
                       S"-XI", S"-XX", S"-XY", S"-XZ",
                       S" YI", S" YX", S" YY", S" YZ",
                       S"-YI", S"-YX", S"-YY", S"-YZ",
                       S" ZI", S" ZX", S" ZY", S" ZZ",
                       S"-ZI", S"-ZX", S"-ZY", S"-ZZ",]                        
        for gate in subtypes(AbstractTwoQubitOperator)
            for stab in stabilizers
                @test apply_inv!(stab, gate(1, 2)) == apply!(stab, inv(CliffordOperator(gate(1, 2), 2)))
            end
        end
    end

    @testset "Apply Inv subsystem order" begin
        state = S"X_Z"
        operation = tCNOT
        indices = [3, 1]

        expected = apply!(copy(state), indices, inv(operation))
        @test apply_inv!(copy(state), indices, operation) == expected
        @test apply_inv!(copy(state), Tuple(indices), operation) == expected

        legacy_state = copy(state)
        @test_deprecated apply_inv!(legacy_state, operation, indices)
        @test legacy_state == expected

        @test_throws MethodError apply!(copy(state), operation, indices)
        @test_throws MethodError apply!(Register(copy(state)), operation, indices)
    end

    @testset "Apply Inv Pauli phases" begin
        for phases in (true, false)
            operation = P"Z"

            state = S"X"
            @test apply_inv!(copy(state), operation; phases) ==
                  apply!(copy(state), operation; phases)

            subsystem_state = S"X_"
            indices = [1]
            @test apply_inv!(copy(subsystem_state), indices, operation; phases) ==
                  apply!(copy(subsystem_state), indices, operation; phases)
        end
    end

end
