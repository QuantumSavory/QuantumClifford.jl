@testitem "Quantum Reed-Muller" begin
    using Test
    using Nemo: echelon_form, matrix, GF
    using QECCore.LinearAlgebra
    using QECCore

    function designed_distance(mat)
        dist = 3
        for row in eachrow(mat)
            count = sum(row)
            if count < dist
                return false
            end
        end
        return true
    end

    @testset "Test QuantumReedMuller(r,m) properties" begin
        for m in 3:10
            H = parity_matrix(QuantumReedMuller(m))
            @test designed_distance(H) == true
            @test code_n(QuantumReedMuller(m)) == 2^m - 1
            @test code_k(QuantumReedMuller(m)) == 1
            @test distance(QuantumReedMuller(m)) == 3
            @test H == parity_matrix(CSS(parity_matrix_x(QuantumReedMuller(m)), parity_matrix_z(QuantumReedMuller(m))))
        end
    end


    @testset "QuantumReedMuller equivalence to simple codes" begin
        # QuantumReedMuller(3) is the Steane7 code.
        # TODO: add function to check if two codes are equivalent, like @test canonicalize!(parity_checks(Steane7())) == parity_checks(QuantumReedMuller(3))
        @test parity_matrix(Steane7()) == parity_matrix(QuantumReedMuller(3))[[3,2,1,6,5,4],:]
        
        # [[15,1,3]] qrm code from table 1 of https://arxiv.org/pdf/1705.0010
        pm = [ 1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
        0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
        0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  1  0  0;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  1  0;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  1  0  0  0  1;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  1  0  1;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  1  1;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  1  0  1;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  1;
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1]
        @test parity_matrix(QuantumReedMuller(4)) == pm
    end
end
