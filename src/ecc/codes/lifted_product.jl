import LinearAlgebra
import Hecke: GroupAlgebra, GroupAlgebraElem, representation_matrix

"""
Lifted product codes [panteleev2021degenerate](@cite) [panteleev2022asymptotically](@cite)

- `A::Matrix{PermGroupRingElem}`: the first base matrix for constructing the lifted product code, whose elements are in a permutation group ring;
- `B::Matrix{PermGroupRingElem}`: the second base matrix for constructing the lifted product code, whose elements are in the same permutation group ring as `A`;
- `repr::Function`: a function that converts the permutation group ring element to a matrix;
 default to be [`permutation_repr`](@ref) for GF(2)-algebra.

A lifted product code is constructed by hypergraph product of the two lifted codes `c₁` and `c₂`.
Here, the hypergraph product is taken over a group ring, which serves as the base ring for both lifted codes.
After the hypergraph product, the parity-check matrices are lifted by `repr`.
The lifting is achieved by applying `repr` to each element of the matrix resulted from the hypergraph product, which is mathematically a linear map from a group algebra element to a binary matrix.

See also: [`LiftedCode`](@ref).
"""
struct LPCode <: AbstractECC
    A::GroupAlgebraElemMatrix
    B::GroupAlgebraElemMatrix
    GA::GroupAlgebra
    repr::Function

    function LPCode(A::GroupAlgebraElemMatrix, B::GroupAlgebraElemMatrix; GA::GroupAlgebra=parent(A[1,1]), repr::Function)
        all(elem.parent == GA for elem in A) && all(elem.parent == GA for elem in B) || error("The base rings of all elements in both matrices must be the same as the group algebra")
        new(A, B, GA, repr)
    end

    function LPCode(c₁::LiftedCode, c₂::LiftedCode; GA::GroupAlgebra=c₁.GA, repr::Function=c₁.repr)
        # we are using the group algebra and the representation function of the first lifted code
        c₁.GA == GA && c₂.GA == GA || error("The base rings of both lifted codes must be the same as the group algebra")
        new(c₁.A, c₂.A, GA, repr)
    end
end

function LPCode(A::FqFieldGroupAlgebraElemMatrix, B::FqFieldGroupAlgebraElemMatrix; GA::GroupAlgebra=parent(A[1,1]))
    LPCode(LiftedCode(A; GA=GA, repr=default_repr), LiftedCode(B; GA=GA, repr=default_repr); GA=GA, repr=default_repr)
end

function LPCode(group_elem_array1::Matrix{<: GroupOrAdditiveGroupElem}, group_elem_array2::Matrix{<: GroupOrAdditiveGroupElem}; GA::GroupAlgebra=group_algebra(GF(2), parent(group_elem_array1[1,1])))
    LPCode(LiftedCode(group_elem_array1; GA=GA), LiftedCode(group_elem_array2; GA=GA); GA=GA, repr=default_repr)
end

function LPCode(shift_array1::Matrix{Int}, shift_array2::Matrix{Int}, l::Int; GA::GroupAlgebra=group_algebra(GF(2), abelian_group(l)))
    LPCode(LiftedCode(shift_array1, l; GA=GA), LiftedCode(shift_array2, l; GA=GA); GA=GA, repr=default_repr)
end

iscss(::Type{LPCode}) = true

function parity_checks_xz(c::LPCode)
    hx, hz = hgp(c.A, c.B')
    hx, hz = lift(c.repr, hx), lift(c.repr, hz)
    return hx, hz
end

parity_checks_x(c::LPCode) = parity_checks_xz(c)[1]

parity_checks_z(c::LPCode) = parity_checks_xz(c)[2]

parity_checks(c::LPCode) = parity_checks(CSS(parity_checks_xz(c)...))

code_n(c::LPCode) = size(c.repr(zero(c.GA)), 2) * (size(c.A, 2) * size(c.B, 1) + size(c.A, 1) * size(c.B, 2))

code_s(c::LPCode) = size(c.repr(zero(c.GA)), 1) * (size(c.A, 1) * size(c.B, 1) + size(c.A, 2) * size(c.B, 2))

"""
Two-block group algebra (2GBA) codes.
"""
function two_block_group_algebra_codes(a::GroupAlgebraElem, b::GroupAlgebraElem)
    A = reshape([a], (1, 1))
    B = reshape([b], (1, 1))
    LPCode(A, B)
end

"""
Generalized bicycle codes.
"""
function generalized_bicycle_codes(a_shifts::Array{Int}, b_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group(l))
    a = sum(GA[n%l+1] for n in a_shifts)
    b = sum(GA[n%l+1] for n in b_shifts)
    two_block_group_algebra_codes(a, b)
end

"""
Bicycle codes.
"""
function bicycle_codes(a_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group(l))
    a = sum(GA[n÷l+1] for n in a_shifts)
    two_block_group_algebra_codes(a, a')
end
