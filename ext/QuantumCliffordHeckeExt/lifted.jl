"""
$TYPEDEF

Classical codes lifted over a group algebra [panteleev2021degenerate](@cite) [panteleev2022asymptotically](@cite).

The parity-check matrix is constructed by applying `repr` to each element of `A`,
which is mathematically a linear map from a group algebra element to a binary matrix.
The size of the parity check matrix will enlarged with each element of `A` being inflated into a matrix.
The procedure is called a lift [panteleev2022asymptotically](@cite).

## Constructors

A lifted code can be constructed via the following approaches:

1. A matrix of group algebra elements.

2. A matrix of group elements, where a group element will be considered as a group algebra element by assigning a unit coefficient.

3. A matrix of integers, where each integer represent the shift of a cyclic permutation. The order of the cyclic permutation should be specified.

The default `GA` is the group algebra of `A[1, 1]`, the default representation is the permutation representation.

## The representation function

In this struct, we use the default representation function `default_repr` to convert a `GF(2)`-group algebra element to a binary matrix.
The default representation, provided by `Hecke`, is the permutation representation.

We also accept a custom representation function.
Such a customization would be useful to reduce the number of bits required by the code construction.

For example, if we use a D4 group for lifting. In our default way, the representation will be `8×8` matrices,
where 8 is the group's order. We can find a `4×4` matrix representation for the group,
with details in [this discussion](https://github.com/QuantumSavory/QuantumClifford.jl/pull/312#discussion_r1682721136).

See also: [`LPCode`](@ref).

$TYPEDFIELDS
"""
struct LiftedCode <: ClassicalCode
    """the base matrix of the code, whose elements are in a group algebra."""
    A::GroupAlgebraElemMatrix
    """the group algebra for which elements in `A` are from."""
    GA::GroupAlgebra
    """
    a function that converts a group algebra element to a binary matrix;
    default to be the permutation representation for GF(2)-algebra."""
    repr::Function

    function LiftedCode(A::GroupAlgebraElemMatrix; GA::GroupAlgebra=parent(A[1, 1]), repr::Function)
        all(elem.parent == GA for elem in A) || error("The base ring of all elements in the code must be the same as the group algebra")
        new(A, GA, repr)
    end
end

default_repr(y::GroupAlgebraElem{FqFieldElem, <: GroupAlgebra}) = Matrix((x -> Bool(Int(lift(ZZ, x)))).(representation_matrix(y)))

"""
The GroupAlgebraElem with `GF(2)` coefficients can be converted to a permutation matrix by `representation_matrix` provided by Hecke.
"""
function LiftedCode(A::Matrix{GroupAlgebraElem{FqFieldElem, <: GroupAlgebra}}; GA::GroupAlgebra=parent(A[1,1]))
    !(characteristic(base_ring(A[1, 1])) == 2) && error("The default permutation representation applies only to GF(2) group algebra; otherwise, a custom representation function should be provided")
    LiftedCode(A; GA=GA, repr=default_repr)
end

function LiftedCode(group_elem_array::Matrix{<: GroupOrAdditiveGroupElem}; GA::GroupAlgebra=group_algebra(GF(2), parent(group_elem_array[1,1])), repr::Union{Function, Nothing}=nothing)
    A = zeros(GA, size(group_elem_array)...)
    for i in axes(group_elem_array, 1), j in axes(group_elem_array, 2)
        A[i, j] = GA[A[i, j]]
    end
    if repr === nothing
        return LiftedCode(A; GA=GA, repr=default_repr)
    else
        return LiftedCode(A; GA=GA, repr=repr)
    end
end

function LiftedCode(shift_array::Matrix{Int}, l::Int; GA::GroupAlgebra=group_algebra(GF(2), abelian_group(l)))
    A = zeros(GA, size(shift_array)...)
    for i in 1:size(shift_array, 1)
        for j in 1:size(shift_array, 2)
            A[i, j] = GA[shift_array[i, j]%l+1]
        end
    end
    return LiftedCode(A; GA=GA, repr=default_repr)
end

function lift(repr::Function, mat::GroupAlgebraElemMatrix)
    vcat([hcat([repr(mat[i, j]) for j in axes(mat, 2)]...) for i in axes(mat, 1)]...)
end

function parity_checks(c::LiftedCode)
    return lift(c.repr, c.A)
end

code_n(c::LiftedCode) = size(c.A, 2) * size(zero(c.GA), 2)

code_s(c::LiftedCode) = size(c.A, 1) * size(zero(c.GA), 1)
