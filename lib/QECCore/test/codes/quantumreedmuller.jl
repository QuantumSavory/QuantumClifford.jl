@testitem "Quantum Reed-Muller" begin
    using Test
    using Nemo: rref, matrix, GF
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
        @test rref(matrix(GF(2),parity_matrix(Steane7()))) == rref(matrix(GF(2),parity_matrix(QuantumReedMuller(3))))
        
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
        @test rref(matrix(GF(2),parity_matrix(QuantumReedMuller(4)))) == rref(matrix(GF(2),pm))
    end
end
