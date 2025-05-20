"""
$TYPEDEF

Lifted product codes ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

A lifted product code is defined by the hypergraph product of a base matrices `A` and the conjugate of another base matrix `B'`.
Here, the hypergraph product is taken over a group algebra, of which the base matrices are consisting.

The binary parity check matrix is obtained by applying `repr` to each element of the matrix resulted from the hypergraph product, which is mathematically a linear map from each group algebra element to a binary matrix.

## Constructors

Multiple constructors are available:

1. Two base matrices of group algebra elements.

2. Two lifted codes, whose base matrices are for quantum code construction.

3. Two base matrices of group elements, where each group element will be considered as a group algebra element by assigning a unit coefficient.

4. Two base matrices of integers, where each integer represent the shift of a cyclic permutation. The order of the cyclic permutation should be specified.

## Examples

A [[882, 24, d ≤ 24]] code from Appendix B of [roffe2023bias](@cite).
We use the 1st constructor to generate the code and check its length and dimension.
During the construction, we do arithmetic operations to get the group algebra elements in base matrices `A` and `B`.
Here `x` is the generator of the group algebra, i.e., offset-1 cyclic permutation, and `GA(1)` is the unit element.

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens; import LinearAlgebra: diagind; using QuantumClifford.ECC;

julia> l = 63; GA = group_algebra(GF(2), abelian_group(l)); x = gens(GA)[];

julia> A = zeros(GA, 7, 7);

julia> A[diagind(A)] .= x^27;

julia> A[diagind(A, -1)] .= x^54;

julia> A[diagind(A, 6)] .= x^54;

julia> A[diagind(A, -2)] .= GA(1);

julia> A[diagind(A, 5)] .= GA(1);

julia> B = reshape([1 + x + x^6], (1, 1));

julia> c1 = LPCode(A, B);

julia> code_n(c1), code_k(c1)
(882, 24)
```

A [[175, 19, d ≤ 10]] code from Eq. (18) in Appendix A of [raveendran2022finite](@cite),
following the 4th constructor.

```jldoctest
julia> import Hecke; using QuantumClifford.ECC;

julia> base_matrix = [0 0 0 0; 0 1 2 5; 0 6 3 1]; l = 7;

julia> c2 = LPCode(base_matrix, l .- base_matrix', l);

julia> code_n(c2), code_k(c2)
(175, 19)
```

## Code subfamilies and convenience constructors for them

- When the base matrices of the `LPCode` are 1×1, the code is called a two-block group-algebra code [`two_block_group_algebra_codes`](@ref).
- When the base matrices of the `LPCode` are 1×1 and their elements are sums of cyclic permutations, the code is called a generalized bicycle code [`generalized_bicycle_codes`](@ref).
- When the two matrices are adjoint to each other, the code is called a bicycle code [`bicycle_codes`](@ref).

## The representation function

We use the default representation function `Hecke.representation_matrix` to convert a `GF(2)`-group algebra element to a binary matrix.
The default representation, provided by `Hecke`, is the permutation representation.

We also accept a custom representation function as detailed in [`LiftedCode`](@ref).

See also: [`LiftedCode`](@ref), [`two_block_group_algebra_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`bicycle_codes`](@ref),
[`haah_cubic_codes`](@ref).

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

# TODO document and doctest example
function LPCode(A::FqFieldGroupAlgebraElemMatrix, B::FqFieldGroupAlgebraElemMatrix; GA::GroupAlgebra=parent(A[1,1]))
    LPCode(LiftedCode(A; GA=GA, repr=representation_matrix), LiftedCode(B; GA=GA, repr=representation_matrix); GA=GA, repr=representation_matrix)
end

# TODO document and doctest example
function LPCode(group_elem_array1::Matrix{<: GroupOrAdditiveGroupElem}, group_elem_array2::Matrix{<: GroupOrAdditiveGroupElem}; GA::GroupAlgebra=group_algebra(GF(2), parent(group_elem_array1[1,1])))
    LPCode(LiftedCode(group_elem_array1; GA=GA), LiftedCode(group_elem_array2; GA=GA); GA=GA, repr=representation_matrix)
end

# TODO document and doctest example
function LPCode(shift_array1::Matrix{Int}, shift_array2::Matrix{Int}, l::Int; GA::GroupAlgebra=group_algebra(GF(2), abelian_group(l)))
    LPCode(LiftedCode(shift_array1, l; GA=GA), LiftedCode(shift_array2, l; GA=GA); GA=GA, repr=representation_matrix)
end

iscss(::Type{LPCode}) = true

function hgp(h₁::GroupAlgebraElemMatrix, h₂::GroupAlgebraElemMatrix)
    r₁, n₁ = size(h₁)
    r₂, n₂ = size(h₂)
    # here we use `permutdims` instead of `transpose` to avoid recursive call
    # convert LinearAlgebra.I to Matrix to fix incompatibility with Julia 1.11.1
    # TODO the performance may be affected by this workaround for large codes
    hx = hcat(kron(h₁, Matrix(LinearAlgebra.I(n₂))), kron(Matrix(LinearAlgebra.I(r₁)), permutedims(group_algebra_conj.(h₂))))
    hz = hcat(kron(Matrix(LinearAlgebra.I(n₁)), h₂), kron(permutedims(group_algebra_conj.(h₁)), Matrix(LinearAlgebra.I(r₂))))
    hx, hz
end

function parity_checks_xz(c::LPCode)
    hx, hz = hgp(c.A, permutedims(group_algebra_conj.(c.B)))
    hx, hz = concat_lift_repr(c.repr,hx), concat_lift_repr(c.repr,hz)
    return hx, hz
end

parity_checks_x(c::LPCode) = parity_checks_xz(c)[1]

parity_checks_z(c::LPCode) = parity_checks_xz(c)[2]

parity_checks(c::LPCode) = parity_checks(CSS(parity_checks_xz(c)...))

code_n(c::LPCode) = size(c.repr(zero(c.GA)), 2) * (size(c.A, 2) * size(c.B, 1) + size(c.A, 1) * size(c.B, 2))

code_s(c::LPCode) = size(c.repr(zero(c.GA)), 1) * (size(c.A, 1) * size(c.B, 1) + size(c.A, 2) * size(c.B, 2))

"""
Two-block group algebra (2BGA) codes, which are a special case of lifted product codes
from two group algebra elements `a` and `b`, used as `1x1` base matrices.

## Examples of 2BGA code subfamilies

### `C₄ x C₂`

Here is an example of a [[56, 28, 2]] 2BGA code from Table 2 of [lin2024quantum](@cite)
with direct product of `C₄ x C₂`.

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens; using QuantumClifford.ECC;

julia> using QuantumClifford.ECC: DistanceMIPAlgorithm; # hide

julia> import HiGHS;

julia> GA = group_algebra(GF(2), abelian_group([14,2]));

julia> x, s = gens(GA);

julia> A = 1 + x^7;

julia> B = 1 + x^7 + s + x^8 + s*x^7 + x;

julia> c = two_block_group_algebra_codes(A,B);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(56, 28, 2)
```

### Bivariate Bicycle codes

Bivariate Bicycle codes are a class of Abelian 2BGA codes formed by the direct product
of two cyclic groups `ℤₗ × ℤₘ`. The parameters `l` and `m` represent the orders of the
first and second cyclic groups, respectively.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/qcga).

A [[756, 16, ≤ 34]] code from Table 3 of [bravyi2024high](@cite):

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens; using QuantumClifford.ECC;

julia> l=21; m=18;

julia> GA = group_algebra(GF(2), abelian_group([l, m]));

julia> x, y = gens(GA);

julia> A = x^3 + y^10 + y^17;

julia> B = y^5 + x^3  + x^19;

julia> c = two_block_group_algebra_codes(A,B);

julia> code_n(c), code_k(c)
(756, 16)
```

### Multivariate Bicycle code

The group algebra of the qubit multivariate bicycle (MB) code with r variables is `𝔽₂[𝐺ᵣ]`,
where `𝐺ᵣ = ℤ/l₁ × ℤ/l₂ × ... × ℤ/lᵣ`.

A [[48, 4, 6]] Weight-6 TB-QLDPC code from Appendix A Table 2 of [voss2024multivariatebicyclecodes](@cite).

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens; using QuantumClifford.ECC;

julia> import QuantumClifford.ECC: DistanceMIPAlgorithm; # hide

julia> import HiGHS;

julia> l=4; m=6;

julia> GA = group_algebra(GF(2), abelian_group([l, m]));

julia> x, y = gens(GA);

julia> z = x*y;

julia> A = x^3 + y^5;

julia> B = x + z^5 + y^5 + y^2;

julia> c = two_block_group_algebra_codes(A, B);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(48, 4, 6)
```

### Coprime Bivariate Bicycle code

The coprime bivariate bicycle (BB) codes are defined by two polynomials `𝑎(𝑥,𝑦)` and `𝑏(𝑥,𝑦)`,
where `𝑙` and `𝑚` are coprime, and can be expressed as univariate polynomials `𝑎(𝜋)` and `𝑏(𝜋)`,
with generator `𝜋 = 𝑥𝑦`. They can be viewed as a special case of Lifted Product construction
based on abelian group `ℤₗ x ℤₘ` where `ℤⱼ` cyclic group of order `j`.

[[108, 12, 6]] coprime-bivariate bicycle (BB) code from Table 2 of [wang2024coprime](@cite).

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens; using QuantumClifford.ECC;

julia> import QuantumClifford.ECC: DistanceMIPAlgorithm; # hide

julia> import HiGHS

julia> l=2; m=27;

julia> GA = group_algebra(GF(2), abelian_group([l*m]));

julia> 𝜋 = gens(GA)[1];

julia> A = 𝜋^2 + 𝜋^5  + 𝜋^44;

julia> B = 𝜋^8 + 𝜋^14 + 𝜋^47;

julia> c = two_block_group_algebra_codes(A, B);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(108, 12, 6)
```

See also: [`LPCode`](@ref), [`generalized_bicycle_codes`](@ref), [`bicycle_codes`](@ref), [`haah_cubic_codes`](@ref).
"""
function two_block_group_algebra_codes(a::GroupAlgebraElem, b::GroupAlgebraElem)
    LPCode([a;;], [b;;])
end

"""
Generalized bicycle codes, which are a special case of *abelian* 2GBA codes (and therefore of lifted product codes).
Here the group is chosen as the cyclic group of order `l`,
and the base matrices `a` and `b` are the sum of the group algebra elements corresponding to the shifts `a_shifts` and `b_shifts`.

See also: [`two_block_group_algebra_codes`](@ref), [`bicycle_codes`](@ref).

A [[254, 28, 14 ≤ d ≤ 20]] code from (A1) in Appendix B of [panteleev2021degenerate](@cite).

```jldoctest
julia> import Hecke; using QuantumClifford.ECC

julia> c = generalized_bicycle_codes([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 127);

julia> code_n(c), code_k(c)
(254, 28)
```

An [[70, 8, 10]] *abelian* 2BGA code from Table 1 of [lin2024quantum](@cite), with cyclic group of
order `l = 35`, illustrates that *abelian* 2BGA codes can be viewed as GB codes.

```jldoctest
julia> import Hecke; using QuantumClifford.ECC

julia> l = 35;

julia> c1 = generalized_bicycle_codes([0, 15, 16, 18], [0, 1, 24, 27], l);

julia> code_n(c1), code_k(c1)
(70, 8)
```
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

See also: [`two_block_group_algebra_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`haah_cubic_codes`](@ref).
""" # TODO doctest example
function bicycle_codes(a_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group(l))
    a = sum(GA[n÷l+1] for n in a_shifts)
    two_block_group_algebra_codes(a, group_algebra_conj(a))
end

"""
Haah’s cubic codes [haah2011local](@cite) can be viewed as generalized bicycle (GB) codes
with the group `G = Cₗ × Cₗ × Cₗ`, where `l` denotes the lattice size. In particular, a GB
code with the group `G = ℤ₃ˣ³` corresponds to a cubic code.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/haah_cubic).

```jldoctest
julia> import Hecke; using QuantumClifford.ECC;

julia> c = haah_cubic_codes([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 6);

julia> code_n(c), code_k(c)
(432, 8)
```

See also: [`bicycle_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`two_block_group_algebra_codes`](@ref).
"""
function haah_cubic_codes(a_shifts::Array{Int}, b_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group([l,l,l]))
    a = sum(GA[n%dim(GA)+1] for n in a_shifts)
    b = sum(GA[n%dim(GA)+1] for n in b_shifts)
    two_block_group_algebra_codes(a, b)
end
