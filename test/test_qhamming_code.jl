using Test
using QuantumClifford
using QuantumClifford.ECC: QHamming

# Define test cases
@testset "QHamming Code Tests" begin
    # Test parity checks generation
    @testset "Parity Checks" begin
        hamming_code = QHamming(3)
        @test size(parity_checks(hamming_code)) == (3, 8)
        @test all(parity_checks(hamming_code) .>= 0)  
        # Add more parity check tests if needed
    end

    # Test block size calculation
    @testset "Block Size" begin
        hamming_code = QHamming(3)
        @test code_n(hamming_code) == 8
        # Add more block size tests if needed
    end
end
