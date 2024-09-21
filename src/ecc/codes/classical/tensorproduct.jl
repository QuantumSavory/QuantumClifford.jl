"""
The classical Tensor-product code was discovered by Wolf in his 1963 paper [wolf1965codes](@cite). It combines two component linear error correcting codes by taking the tensor product of their parity check matrices.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/tensor).
"""
struct TensorProduct
    A::AbstractMatrix
    B::AbstractMatrix
    function TensorProduct(A, B)
        all(x -> x isa Bool, A) && all(x -> x isa Bool, B) || throw(ArgumentError("All matrices should be defined over the Boolean field, GF(2)"))
        new(A, B)
    end
end

function parity_checks(tp::TensorProduct)
    A, B = tp.A, tp.B
    C = kronecker_product(matrix(ZZ, B), matrix(ZZ, A))
    C = Matrix{Bool}(Matrix{Int}(C))
    return C
end
