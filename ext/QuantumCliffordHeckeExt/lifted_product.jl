"""
$TYPEDEF

Lifted product codes [panteleev2021degenerate](@cite) [panteleev2022asymptotically](@cite)

A lifted product code is defined by the hypergraph product of a base matrices `A` and the conjugate of another base matrix `B'`.
Here, the hypergraph product is taken over a group algebra, of which the base matrices are consisting.

The parity check matrix is obtained by applying `repr` to each element of the matrix resulted from the hypergraph product, which is mathematically a linear map from each group algebra element to a binary matrix.

## Constructors

1. Two base matrices of group algebra elements.

2. Two lifted codes, whose base matrices are for quantum code construction. 

3. Two base matrices of group elements, where each group element will be considered as a group algebra element by assigning a unit coefficient.

4. Two base matrices of integers, where each integer represent the shift of a cyclic permutation. The order of the cyclic permutation should be specified.

## The representation function

In this struct, we use the default representation function `default_repr` to convert a `GF(2)`-group algebra element to a binary matrix.
The default representation, provided by `Hecke`, is the permutation representation.

We also accept a custom representation function. The reasons are detailed in [`LiftedCode`](@ref).

See also: [`LiftedCode`](@ref), [`two_block_group_algebra_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`bicycle_codes`](@ref).

$TYPEDFIELDS
"""
struct LPCode <: AbstractECC
    """the first base matrix of the code, whose elements are in a group algebra."""
    A::GroupAlgebraElemMatrix
    """the second base matrix of the code, whose elements are in the same group algebra as `A`."""
    B::GroupAlgebraElemMatrix
    """the group algebra for which elements in `A` and `B` are from."""
    GA::GroupAlgebra
    """
    a function that converts a group algebra element to a binary matrix;
    default to be the permutation representation for GF(2)-algebra."""
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
Two-block group algebra (2GBA) codes, which are a special case of lifted product codes
from two group algebra elements `a` and `b`, used as `1x1` base matrices.

See also: [`LPCode`](@ref), [`generalized_bicycle_codes`](@ref), [`bicycle_codes`](@ref)
"""
function two_block_group_algebra_codes(a::GroupAlgebraElem, b::GroupAlgebraElem)
    A = reshape([a], (1, 1))
    B = reshape([b], (1, 1))
    LPCode(A, B)
end

"""
Generalized bicycle codes, which are a special case of 2GBA codes (and therefore of lifted product codes).
Here the group is choosen as the cyclic group of order `l`,
and the base matrices `a` and `b` are the sum of the group algebra elements corresponding to the shifts `a_shifts` and `b_shifts`.

See also: [`two_block_group_algebra_codes`](@ref), [`bicycle_codes`](@ref).
"""
function generalized_bicycle_codes(a_shifts::Array{Int}, b_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group(l))
    a = sum(GA[n%l+1] for n in a_shifts)
    b = sum(GA[n%l+1] for n in b_shifts)
    two_block_group_algebra_codes(a, b)
end

"""
Bicycle codes are a special case of generalized bicycle codes,
where `a` and `b` are conjugate to each other.
The order of the cyclic group is `l`, and the shifts `a_shifts` and `b_shifts` are reverse to each other.

See also: [`two_block_group_algebra_codes`](@ref), [`generalized_bicycle_codes`](@ref).
"""
function bicycle_codes(a_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group(l))
    a = sum(GA[n÷l+1] for n in a_shifts)
    two_block_group_algebra_codes(a, a')
end
