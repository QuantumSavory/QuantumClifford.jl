"""Hypergraph product of two classical codes."""
function hgp(h₁,h₂)
    r₁, n₁ = size(h₁)
    r₂, n₂ = size(h₂)
    hx = hcat(kron(h₁, LinearAlgebra.I(n₂)), kron(LinearAlgebra.I(r₁), h₂'))
    hz = hcat(kron(LinearAlgebra.I(n₁), h₂), kron(h₁', LinearAlgebra.I(r₂)))
    hx, hz
end

import AbstractAlgebra: Group, GroupElem, AdditiveGroup, AdditiveGroupElem
import Hecke: GroupAlgebra, GroupAlgebraElem, dim, base_ring, multiplication_table, coefficients
import Base: adjoint
import LinearAlgebra

"""
Compute the adjoint of a group algebra element.
The adjoint is defined as the conjugate of the element in the group algebra,
i.e. the inverse of the element in the associated group.
"""
function adjoint(a::GroupAlgebraElem{T}) where T
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

"""
The difference between Group and AdditiveGroup is that the former one uses * operations while the latter uses + operations;
they are all supported by Hecke.GroupAlgebraElem
"""
# const GroupOrAdditiveGroup = Union{Group,AdditiveGroup}

const GroupOrAdditiveGroupElem = Union{GroupElem,AdditiveGroupElem}

const GroupAlgebraElemMatrix = Union{
    Matrix{<:GroupAlgebraElem},
    LinearAlgebra.Adjoint{<:GroupAlgebraElem,<:Matrix{<:GroupAlgebraElem}}
}

const FqFieldGroupAlgebraElemMatrix = Union{
    Matrix{<:GroupAlgebraElem{FqFieldElem,<:GroupAlgebra}},
    LinearAlgebra.Adjoint{<:GroupAlgebraElem{FqFieldElem,<:GroupAlgebra},<:Matrix{<:GroupAlgebraElem{FqFieldElem,<:GroupAlgebra}}}
}
