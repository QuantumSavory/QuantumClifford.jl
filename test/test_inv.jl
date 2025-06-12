@testitem "Apply Inv" begin
    using InteractiveUtils

    @testset "Apply Inv" begin
        stabilizers = [S"I", S"X", S"Y", S"Z"]
        for gate in subtypes(AbstractSingleQubitOperator)
            gate == SingleQubitOperator && continue
            for stab in stabilizers
                @test apply_inv!(stab, gate(1)) == apply!(stab, inv(CliffordOperator(gate(1), 1)))
            end
        end
    end

    @testset "Apply Inv two-qubit" begin
        stabilizers = [S"II", S"IX", S"IY", S"IZ",
                       S"XI", S"XX", S"XY", S"XZ",
                       S"YI", S"YX", S"YY", S"YZ",
                       S"ZI", S"ZX", S"ZY", S"ZZ"]
        for gate in subtypes(AbstractTwoQubitOperator)
            for stab in stabilizers
                @test apply_inv!(stab, gate(1, 2)) == apply!(stab, inv(CliffordOperator(gate(1, 2), 2)))
            end
        end
    end
    
end