"""
$TYPEDEF

Classical codes lifted over a group algebra, used for lifted product code construction ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

The parity-check matrix is constructed by applying `repr` to each element of `A`,
which is mathematically a linear map from a group algebra element to a binary matrix.
The size of the parity check matrix will enlarged with each element of `A` being inflated into a matrix.
The procedure is called a lift [panteleev2022asymptotically](@cite).

## Constructors

A lifted code can be constructed via the following approaches:

1. A matrix of group algebra elements.

2. A matrix of group elements, where a group element will be considered as a group algebra element by assigning a unit coefficient.

3. A matrix of integers, where each integer represent the shift of a cyclic permutation. The order of the cyclic permutation should be specified.

The default `GA` is the group algebra of `A[1, 1]`, the default representation `repr` is the permutation representation.

Below is a list of all constructors:

$METHODLIST

## The representation function `repr`

We use the default representation function `Hecke.representation_matrix` to convert a `GF(2)`-group algebra element to a binary matrix.
The default representation, provided by `Hecke`, is the permutation representation.

We also accept a custom representation function (the `repr` field of the constructor).
Whatever the representation, the matrix elements need to be convertible to Integers (e.g. permit `lift(ZZ, ...)`).
Such a customization would be useful to reduce the number of bits required by the code construction.

For example, if we use a D4 group for lifting, our default representation will be `8×8` permutation matrices,
where 8 is the group's order.
However, we can find a `4×4` matrix representation for the group,
e.g. by using the typical [`2×2` representation](https://en.wikipedia.org/wiki/Dihedral_group)
and converting it into binary representation by replacing "1" with the Pauli I, and "-1" with the Pauli X matrix.

See also: [`QuantumClifford.ECC.LPCode`](@ref).

All fields:

$TYPEDFIELDS
"""
struct LiftedCode <: AbstractCECC
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

#"""
#`LiftedCode` constructor using the default `GF(2)` representation (coefficients converted to a permutation matrix by `representation_matrix` provided by Hecke).
#"""
# TODO doctest example and document the other constructors below
function LiftedCode(A::Matrix{GroupAlgebraElem{FqFieldElem, <: GroupAlgebra}}; GA::GroupAlgebra=parent(A[1,1]))
    !(characteristic(base_ring(A[1, 1])) == 2) && error("The default permutation representation applies only to GF(2) group algebra; otherwise, a custom representation function should be provided")
    LiftedCode(A; GA=GA, repr=representation_matrix)
end

"""
Constructs a group algebra code by embedding a matrix of group elements into the
specified group algebra `GA`, with optional custom representation  `repr`.

# Example

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens, representation_matrix

julia> import QuantumClifford.ECC: LiftedCode, code_n, code_k, code_s

julia> import QuantumClifford.ECC.QECCore: parity_matrix

julia> l = 12; GA = group_algebra(GF(2), abelian_group(l)); x = gens(GA)[];

julia> B = reshape([1 + x + x^3 + x^6], (1, 1));

julia> c = LiftedCode(B, repr = representation_matrix);

julia> code = parity_matrix(c)
12×12 Matrix{Bool}:
 1  0  0  0  0  0  1  0  0  1  0  1
 1  1  0  0  0  0  0  1  0  0  1  0
 0  1  1  0  0  0  0  0  1  0  0  1
 1  0  1  1  0  0  0  0  0  1  0  0
 0  1  0  1  1  0  0  0  0  0  1  0
 0  0  1  0  1  1  0  0  0  0  0  1
 1  0  0  1  0  1  1  0  0  0  0  0
 0  1  0  0  1  0  1  1  0  0  0  0
 0  0  1  0  0  1  0  1  1  0  0  0
 0  0  0  1  0  0  1  0  1  1  0  0
 0  0  0  0  1  0  0  1  0  1  1  0
 0  0  0  0  0  1  0  0  1  0  1  1

julia> code_n(c), code_k(c), code_s(c)
(12, 3, 12)
```
"""
function LiftedCode(group_elem_array::Matrix{<: GroupOrAdditiveGroupElem}; GA::GroupAlgebra=group_algebra(GF(2), parent(group_elem_array[1,1])), repr::Union{Function, Nothing}=nothing)
    A = zeros(GA, size(group_elem_array)...)
    for i in axes(group_elem_array, 1), j in axes(group_elem_array, 2)
        A[i, j] = GA[A[i, j]]
    end
    if repr === nothing
        return LiftedCode(A; GA=GA, repr=representation_matrix)
    else
        return LiftedCode(A; GA=GA, repr=repr)
    end
end

"""
Constructs a group algebra code over `GF(2)` by lifting a matrix of cyclic shifts
(entries modulo `l`)  to the group algebra of the abelian group `ℤ/lℤ` of order `l`.

# Example

```jldoctest
julia> import Hecke; import QuantumClifford.ECC: LiftedCode, code_n, code_k, code_s

julia> import QuantumClifford.ECC.QECCore: parity_matrix

julia> base_matrix = [0 0 0 0; 0 1 2 5; 0 6 3 1]; l = 3;

julia> c = LiftedCode(base_matrix, l);

julia> code = parity_matrix(c)
9×12 Matrix{Bool}:
 1  0  0  1  0  0  1  0  0  1  0  0
 0  1  0  0  1  0  0  1  0  0  1  0
 0  0  1  0  0  1  0  0  1  0  0  1
 1  0  0  0  0  1  0  1  0  0  1  0
 0  1  0  1  0  0  0  0  1  0  0  1
 0  0  1  0  1  0  1  0  0  1  0  0
 1  0  0  1  0  0  1  0  0  0  0  1
 0  1  0  0  1  0  0  1  0  1  0  0
 0  0  1  0  0  1  0  0  1  0  1  0

julia> code_n(c), code_k(c), code_s(c)
(12, 5, 9)
```
"""
function LiftedCode(shift_array::Matrix{Int}, l::Int; GA::GroupAlgebra=group_algebra(GF(2), abelian_group(l)), repr=representation_matrix)
    A = zeros(GA, size(shift_array)...)
    for i in 1:size(shift_array, 1)
        for j in 1:size(shift_array, 2)
            A[i, j] = GA[shift_array[i, j]%l+1]
        end
    end
    return LiftedCode(A; GA, repr)
end

lift_to_bool(x) = Bool(Int(lift(ZZ,x)))

function concat_lift_repr(repr, mat)
    x = repr.(mat)
    y = hvcat(size(x,2), transpose(x)...)
    z = Matrix(lift_to_bool.(y))
    return z
end

function parity_matrix(c::LiftedCode)
    return concat_lift_repr(c.repr, c.A)
end

code_n(c::LiftedCode) = size(c.A, 2) * order(group(c.GA))

code_s(c::LiftedCode) = size(c.A, 1) * order(group(c.GA))

code_k(c::LiftedCode) = code_n(c) - rank(matrix(GF(2), parity_matrix(c)))
