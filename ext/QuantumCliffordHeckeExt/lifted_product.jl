right_repr_matrix(x) = representation_matrix(x, :right)
left_repr_matrix(x) = representation_matrix(x, :left)

"""
$TYPEDEF

Lifted product codes ([panteleev2021degenerate](@cite), [panteleev2022asymptotically](@cite))

A lifted product code is defined by the hypergraph product of a base matrices `A` and the conjugate of another base matrix `B'`.
Here, the hypergraph product is taken over a group algebra, of which the base matrices are consisting.

The binary parity check matrices are obtained by applying `A_repr` and `B_repr` representation maps to each element of the base matrices. These linear transformations convert group algebra elements to their matrix representations while preserving the CSS orthogonality condition.

## Mathematical Framework

Given classical parity-check matrices:

- ``A \\in \\mathbb{F}_q^{m_a \\times n_a}``

- ``B \\in \\mathbb{F}_q^{m_b \\times n_b}``

The lifted product construction produces quantum CSS codes with parity-check matrices:

```math
\\begin{aligned}
    H_X &= [A \\otimes I_{m_b}, -I_{m_a} \\otimes B] \\\\
    H_Z &= [I_{n_a} \\otimes B^*, A^* \\otimes I_{n_b}]
\\end{aligned}

### Commutative Algebra

When `R` is *commutative*, a single representation suffices since all elements naturally commute. Here ``\\rho(a) = \\lambda(a)`` for all ``a \\in R``.

### Non-Commutative Algebra

When `R` is *non-commutative*, distinct representations are essential:

- `A_repr` implements the right regular representation: ``\\rho(a)x = xa``

- `B_repr` implements the left regular representation: ``\\lambda(b)x = bx``

These ensure the critical commutation relation:

```math
\\begin{aligned}
    \\rho(a)\\lambda(b) = \\lambda(b)\\rho(a)
\\end{aligned}
```

which follows from the associative property:

```math
\\begin{aligned}
    \\rho(a)\\lambda(b)(x) = b(xa) = (bx)a = \\lambda(b)\\rho(a)(x)
\\end{aligned}
```

## Constructors

Multiple constructors are available:

1. Two base matrices of group algebra elements.

2. Two lifted codes, whose base matrices are for quantum code construction.

3. Two base matrices of group elements, where each group element will be considered as a group algebra element by assigning a unit coefficient.

4. Two base matrices of integers, where each integer represent the shift of a cyclic permutation. The order of the cyclic permutation should be specified.

Below is a list of all constructors:

$METHODLIST

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

All fields:

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
    a function that converts a group algebra element to a binary matrix for `A`;
    default to be the right regular representation for `GF(2)`-algebra."""
    A_repr::Function
    """
    a function that converts a group algebra element to a binary matrix for B;
    default to be the left regular representation for `GF(2)`-algebra."""
    B_repr::Function

    # TODO document and doctest example
    function LPCode(
        A::GroupAlgebraElemMatrix,
        B::GroupAlgebraElemMatrix;
        GA::GroupAlgebra=parent(A[1,1]),
        A_repr::Function=right_repr_matrix,
        B_repr::Function=left_repr_matrix,
        repr::Union{Function,Nothing}=nothing
    )
        if repr !== nothing # override the default A_repr/B_repr (exists for backward compat)
            is_commutative(GA) || throw(ArgumentError("The group algebra must be commutative when using a single `repr` function, which is not the case here. Please specify separate `A_repr` and `B_repr` instead of a single `repr`. The default choice of `A_repr=right_repr_matrix, B_repr=left_repr_matrix` is frequently sufficient."))
            A_repr = B_repr = repr
        end
        all(elem.parent == GA for elem in A) && all(elem.parent == GA for elem in B) || error("The base rings of all elements in both matrices must be the same as the group algebra")
        new(A, B, GA, A_repr, B_repr)
    end

    # TODO document and doctest example
    function LPCode(
        c₁::LiftedCode,
        c₂::LiftedCode;
        GA::GroupAlgebra=c₁.GA,
        A_repr::Function=c₁.repr,
        B_repr::Function=c₂.repr,
        repr::Union{Function,Nothing}=nothing
    )
        if repr !== nothing # override the default A_repr/B_repr (exists for backward compat)
            is_commutative(GA) || throw(ArgumentError("The group algebra must be commutative when using a single `repr` function, which is not the case here. Please specify separate `A_repr` and `B_repr` instead of a single `repr`. The default choice of `A_repr=right_repr_matrix, B_repr=left_repr_matrix` is frequently sufficient."))
            A_repr = B_repr = repr
        end
        # We are using the representation function of each lifted code.
        # We are using their group algebras as well (asserting that they are the same).
        c₁.GA == GA && c₂.GA == GA || error("The base rings of both lifted codes must be the same as the group algebra")
        new(c₁.A, c₂.A, GA, A_repr, B_repr)
    end
end

# TODO document and doctest example
function LPCode(A::FqFieldGroupAlgebraElemMatrix, B::FqFieldGroupAlgebraElemMatrix; GA::GroupAlgebra=parent(A[1,1]))
    LPCode(
        LiftedCode(A; GA=GA, repr=right_repr_matrix),
        LiftedCode(B; GA=GA, repr=left_repr_matrix)
    )
end

# TODO document and doctest example
function LPCode(group_elem_array1::Matrix{<: GroupOrAdditiveGroupElem}, group_elem_array2::Matrix{<: GroupOrAdditiveGroupElem}; GA::GroupAlgebra=group_algebra(GF(2), parent(group_elem_array1[1,1])))
    LPCode(
        LiftedCode(group_elem_array1; GA=GA, repr=right_repr_matrix),
        LiftedCode(group_elem_array2; GA=GA, repr=left_repr_matrix)
    )
end

# TODO document and doctest example
function LPCode(shift_array1::Matrix{Int}, shift_array2::Matrix{Int}, l::Int; GA::GroupAlgebra=group_algebra(GF(2), abelian_group(l)))
    LPCode(
        LiftedCode(shift_array1, l; GA=GA, repr=right_repr_matrix),
        LiftedCode(shift_array2, l; GA=GA, repr=left_repr_matrix)
    )
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

function parity_matrix_xz(c::LPCode)
    A, B = c.A, c.B
    ma, na = size(A)
    mb, nb = size(B)
    # hx = [A ⊗ I  | I  ⊗ B]
    # hz = [I ⊗ B* | A* ⊗ I]
    hx_raw, hz_raw = hgp(A, permutedims(group_algebra_conj.(B)))
    hx_b₁cols = hz_b₁cols = na*mb # size(A ⊗ I, 2) == size(I ⊗ B*, 2)
    hx = hcat(
        concat_lift_repr(c.A_repr, hx_raw[:, 1:hx_b₁cols]), # ρ(A ⊗ I)
        concat_lift_repr(c.B_repr, hx_raw[:, hx_b₁cols+1:end]) # λ(I ⊗ B)
    )
    hz = hcat(
        concat_lift_repr(c.B_repr, hz_raw[:, 1:hz_b₁cols]), # λ(I ⊗ B*)
        concat_lift_repr(c.A_repr, hz_raw[:, hz_b₁cols+1:end]) # ρ(A* ⊗ I)
    )
    # hx = [ρ(A ⊗  I) | λ(I  ⊗ B)]
    # hz = [λ(I ⊗ B*) | ρ(A* ⊗ I)]
    return hx, hz
end

"""A simpler (older) implementation of the parity check matrix representation routine, used purely for testing. Valid only for commuting algebras."""
function _parity_matrix_xz_if_comm_algebra(c::LPCode)
    hx, hz = hgp(c.A, permutedims(group_algebra_conj.(c.B)))
    hx, hz = concat_lift_repr(c.A_repr,hx), concat_lift_repr(c.A_repr,hz)
    return hx, hz
end

parity_matrix_x(c::LPCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::LPCode) = parity_matrix_xz(c)[2]

parity_checks(c::LPCode) = parity_checks(CSS(parity_matrix_xz(c)...))

code_n(c::LPCode) = size(c.A_repr(zero(c.GA)), 2) * (size(c.A, 2) * size(c.B, 1) + size(c.A, 1) * size(c.B, 2))

code_s(c::LPCode) = size(c.A_repr(zero(c.GA)), 1) * (size(c.A, 1) * size(c.B, 1) + size(c.A, 2) * size(c.B, 2))

"""
$TYPEDSIGNATURES

Two-block group algebra (2BGA) codes, which are a special case of lifted product codes
from two group algebra elements `a` and `b`, used as `1×1` base matrices.
To build them, you pick a group and specific generators for that group,
then you pick two polynomials made of the group generators,
and then, behind the scenes, these two polynomials `a` and `b` are piped
to the lifted product code constructor as the elements of `1×1` matrices.

See also: [`LPCode`](@ref), [`generalized_bicycle_codes`](@ref), [`bicycle_codes`](@ref), [`haah_cubic_codes`](@ref),
[`honeycomb_color_codes`](@ref).

## Examples of 2BGA code subfamilies

### `C₄ x C₂`

Here is an example of a [[56, 28, 2]] 2BGA code from Table 2 of [lin2024quantum](@cite)
build out of polymonials of generators of the direct product `C₄ × C₂`.

```jldoctest
julia> import Hecke: group_algebra, GF, abelian_group, gens; using QuantumClifford.ECC;

julia> GA = group_algebra(GF(2), abelian_group([14,2]));

julia> x, s = gens(GA);

julia> A = 1 + x^7;

julia> B = 1 + x^7 + s + x^8 + s*x^7 + x;

julia> c = two_block_group_algebra_codes(A,B);

julia> import HiGHS;

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

julia> l=4; m=6;

julia> GA = group_algebra(GF(2), abelian_group([l, m]));

julia> x, y = gens(GA);

julia> z = x*y;

julia> A = x^3 + y^5;

julia> B = x + z^5 + y^5 + y^2;

julia> c = two_block_group_algebra_codes(A, B);

julia> import HiGHS

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

julia> l=2; m=27;

julia> GA = group_algebra(GF(2), abelian_group([l*m]));

julia> 𝜋 = gens(GA)[1];

julia> A = 𝜋^2 + 𝜋^5  + 𝜋^44;

julia> B = 𝜋^8 + 𝜋^14 + 𝜋^47;

julia> c = two_block_group_algebra_codes(A, B);

julia> import HiGHS

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(108, 12, 6)
```
"""
function two_block_group_algebra_codes(a::GroupAlgebraElem, b::GroupAlgebraElem)
    LPCode([a;;], [b;;])
end

"""
$TYPEDSIGNATURES

Generalized bicycle codes, which are a special case of *abelian* 2GBA codes (and therefore of lifted product codes).
Here the group is chosen as the cyclic group of order `l`,
and the base matrices `a` and `b` are the sum of the group algebra elements corresponding to the shifts `a_shifts` and `b_shifts`.

Behind the scenes, the shifts are converted to the corresponding group algebra elements and piped to [`two_block_group_algebra_codes`](@ref).

See also: [`two_block_group_algebra_codes`](@ref), [`bicycle_codes`](@ref).

## Examples

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
$TYPEDSIGNATURES

Bicycle codes are a special case of generalized bicycle codes,
where `a` and `b` are conjugate to each other.
The order of the cyclic group is `l`, and the shifts `a_shifts` and `b_shifts` are reverse to each other.
Thus you need to provide only the `a_shifts` and the rest of the conversions and conjugations are taken care of.

See also: [`two_block_group_algebra_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`haah_cubic_codes`](@ref).
""" # TODO doctest example
function bicycle_codes(a_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group(l))
    a = sum(GA[n÷l+1] for n in a_shifts)
    two_block_group_algebra_codes(a, group_algebra_conj(a))
end

"""
$TYPEDSIGNATURES

Haah’s cubic codes [haah2011local](@cite) can be viewed as generalized bicycle (GB) codes
with the group `G = Cₗ × Cₗ × Cₗ`, where `l` denotes the lattice size. In particular, a GB
code with the group `G = ℤ₃ˣ³` corresponds to a cubic code.

Behind the scenes, this function is just a simple shortcut for preparing the group `G`,
before piping the arguments to [`generalized_bicycle_codes`](@ref).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/haah_cubic).

See also: [`bicycle_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`two_block_group_algebra_codes`](@ref).

## Examples

```jldoctest
julia> import Hecke; using QuantumClifford.ECC;

julia> c = haah_cubic_codes([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 6);

julia> code_n(c), code_k(c)
(432, 8)
```
"""
function haah_cubic_codes(a_shifts::Array{Int}, b_shifts::Array{Int}, l::Int)
    GA = group_algebra(GF(2), abelian_group([l,l,l]))
    a = sum(GA[n%dim(GA)+1] for n in a_shifts)
    b = sum(GA[n%dim(GA)+1] for n in b_shifts)
    two_block_group_algebra_codes(a, b)
end

"""
Haah’s cubic code is defined as ``\\text{LP}(1 + x + y + z, 1 + xy + xz + yz)``
where ``\\text{LP}`` is the lifted product code, and `x`, `y`, `z` are elements
of the ring ``R = \\mathbb{F}_2[x, y, z] / (x^L - 1, y^L - 1, z^L - 1)``. Here
``\\mathbb{F}_2`` is the finite field of order `2` and `L` is the lattice size.
The ring ``R`` is the group algebra ``\\mathbb{F}_qG`` of a finite group `G`, where
``G = (C_L)^3`` and ``C_L`` is the cyclic group of order `L`. This method of Haah's
cubic code construction is outlined in Appendix B of [panteleev2022asymptotically](@cite).

Here is an example of a `[[1024, 30, 13 ≤ d ≤ 32]]` Haah's cubic code from Appendix B,
code D of [panteleev2021degenerate](@cite) on the `8 × 8 × 8` Lattice.

```jldoctest
julia> import Hecke; using QuantumClifford.ECC;

julia> l = 8;

julia> c = haah_cubic_codes(l);

julia> code_n(c), code_k(c)
(1024, 30)
```

See also: [`bicycle_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`two_block_group_algebra_codes`](@ref),
[`honeycomb_color_codes`](@ref).
"""
function haah_cubic_codes(l::Int)
    GA = group_algebra(GF(2), abelian_group([l,l,l]))
    x, y, z = gens(GA)
    c = [1 + x + y + z;;]
    d = [1 + x*y + x*z + y*z;;]
    LPCode(c,d)
end

"""
The honeycomb color codes [eberhardt2024logical](@cite) are exactly the Bivariate
Bicycle (BB) codes defined by the polynomials `c = 1 + x + xy` and `d = 1 + y + xy`,
provided that both `ℓ` and `m` are divisible by three. This `6.6.6` code is an example of BB
code, as it represents a special case.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/triangular_color).

```jldoctest
julia> import Hecke; using QuantumClifford.ECC;

julia> ℓ = 9; m = 6;

julia> c = honeycomb_color_codes(ℓ, m);

julia> code_n(c), code_k(c)
(108, 4)
```

See also: [`bicycle_codes`](@ref), [`generalized_bicycle_codes`](@ref), [`two_block_group_algebra_codes`](@ref),
[`honeycomb_color_codes`](@ref).
"""
function honeycomb_color_codes(ℓ::Int, m::Int)
    (ℓ % 3 == 0 && m % 3 == 0) || throw(ArgumentError("Both ℓ and m must be divisible by 3"))
    GA = group_algebra(GF(2), abelian_group([ℓ, m]))
    x, y = gens(GA)
    c = 1 + x + x*y
    d = 1 + y + x*y
    two_block_group_algebra_codes(c, d)
end
