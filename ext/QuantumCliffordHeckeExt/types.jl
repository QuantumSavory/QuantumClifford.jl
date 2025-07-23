const GroupOrAdditiveGroupElem = Union{GroupElem,AdditiveGroupElem}

const GroupAlgebraElemMatrix = Union{
    Matrix{<:GroupAlgebraElem},
    LinearAlgebra.Adjoint{<:GroupAlgebraElem,<:Matrix{<:GroupAlgebraElem}}
}

const FqFieldGroupAlgebraElemMatrix = Union{
    Matrix{<:GroupAlgebraElem{FqFieldElem,<:GroupAlgebra}},
    LinearAlgebra.Adjoint{<:GroupAlgebraElem{FqFieldElem,<:GroupAlgebra},<:Matrix{<:GroupAlgebraElem{FqFieldElem,<:GroupAlgebra}}}
}

"""
Compute the conjugate of a group algebra element.
The conjugate is defined by inversing elements in the associated group.
"""
function group_algebra_conj(a::GroupAlgebraElem{T}) where T
    A = parent(a)
    d = dim(A)
    v = Vector{T}(undef, d)
    for i in 1:d
      v[i] = zero(base_ring(A))
    end
    id_index = findfirst(x -> x == 1, one(A).coeffs)
    # t = zero(base_ring(A))
    mt = multiplication_table(A, copy = false)
    acoeff = coefficients(a, copy = false)
    for i in 1:d
        if acoeff[i] != 0
            k = findfirst(x -> x==id_index, mt[i, :]) # find the inverse of i-th element in the group
            v[k] += acoeff[i]
        end
    end
    return A(v)
end
