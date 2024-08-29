@testitem "Tensor Producr codes" begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: TensorProduct, parity_checks
    A = [1 1 1 1 1] .== 1 
    B = [1 1 1 0 1 0 0; 1 1 0 1 0 1 0; 1 0 1 1 0 0 1] .== 1
    # Example 1 of [wolf2006introduction](@cite)
    @test parity_checks(TensorProduct(A, B)) == Matrix{Bool}([ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0;
                                                               1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0;
                                                               1 1 1 1 1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1])
end
