"""
The binary quantum Tensor Product Code(TQC) is derived from classical TPCs using the CSS construction. A classical TPC is dual-containing if at least one of its component codes is.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_tensor_product).
"""
struct QuantumTensorProduct
    A::AbstractMatrix
    B::AbstractMatrix
    function QuantumTensorProduct(A, B)
        all(x -> x isa Bool, A) && all(x -> x isa Bool, B) || throw(ArgumentError("All matrices should be defined over the Boolean field, GF(2)."))
        new(A, B)
    end
end

function parity_checks(Tₚ::QuantumTensorProduct)
    T = parity_checks(TensorProduct(Tₚ.A, Tₚ.B))
    H = Stabilizer(CSS(T, T))
    return H
end
