"""
    $TYPEDEF

A generalization of the circulant bivariate bicycle quantum LDPC code that allows
for more than `3` terms in the ``A`` and `B`` matrices. The original construction was
introduced in [bravyi2024high](@cite) specifically employed three-term polynomials
where ``A = A_1 + A_2 + A_3`` and ``B = B_1 + B_2 + B_3``, with each ``A_i`` and ``B_i``
representing either powers of ``x = S_ℓ ⊗ I_m`` or ``y = I_ℓ ⊗ S_m``.

The construction uses identity and cyclic shift matrices where `Iₗ` is the `l × l`
identity matrix  and `Sₗ` is the cyclic shift matrix of the same size. We use the
matrices `x = Sₗ ⊗ Iₘ` and  `y = Iₗ ⊗ Sₘ`. The generalized code is represented by
matrices `A` and `B`, defined as sums of terms where each term is a power of `x` or
`y`. The check matrices are: `Hx = [A | B]` and  `Hz = [B'|A']`.

## Circulant Bivariate Bicycle Code

The bivariate bicycle (BB) quantum LDPC code was introduced in [bravyi2024high](@cite).
This construction uses identity and cyclic shift matrices because we consider `Iₗ`
as the `l × l` identity matrix and `Sₗ` as the cyclic shift matrix of the same size,
where each row of `Sₗ` has a single '1' at the column `(i + 1) mod l`. We use the
matrices `x = Sₗ ⊗ Iₘ` and `y = Iₗ ⊗ Sₘ`. The BB code is represented by matrices
`A` and `B`, defined as: `A = A₁ + A₂ + A₃` and `B = B₁ + B₂ + B₃`. The addition and
multiplication operations on binary matrices are performed modulo 2. The check matrices
are: `Hx = [A|B]` and `Hz = [B'|A']`. Both `Hx` and `Hz` are `(n/2)×n` matrices.

!!! note
    Both A and B are matrices in which each row and column contains exactly three non-zero
    entries when using a 3-term polynomial representation [bravyi2024high](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/qcga).

### Example

Here is `[[756, 16, ≤ 34]]` circulant bivariate bicycle code from Table 1 of [bravyi2024high](@cite)
with polynomials ``A = x^3 + y^10 + y^17`` and ``B = y^5 + x^3 + x^19``. 

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> l, m = 21, 18;

julia> A = [(:x,3), (:y,10), (:y,17)];

julia> B = [(:y,5), (:x,3), (:x,19)];

julia> c = GeneralizedCirculantBivariateBicycle(l, m, A, B);

julia> code_n(c), code_k(c)
(756, 16)
```

## Generalized Circulant Bivariate Bicycle Code

Here is an example of  `[[128, 14, 12]]` generalized circulant bivariate bicycle code
that uses 4-term polynomials ``c = x^2 + y + y^3 + y^4`` and ``d = y^2 + x + x^3 + x^4``
with group orders ``l, m = 8, 8`` from [eberhardt2024logical](@cite).

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> l, m = 8, 8;

julia> A = [(:x,2), (:y,1), (:y,3), (:y,4)];

julia> B = [(:y,2), (:x,1), (:x,3), (:x,4)];

julia> c = GeneralizedCirculantBivariateBicycle(l, m, A, B);

julia> code_n(c), code_k(c)
(128, 14)
```

!!! note
    The Bivariate Bicycle code ``\\mathrm{QC}(A,B)`` is a specific instance of the
    *Lifted Product* construction, where the underlying group is the direct product 
    ``\\mathbb{Z}_\\ell \\times \\mathbb{Z}_m`` (with ``\\mathbb{Z}_j`` denoting the
    cyclic group of order ``j``).

### Fields
    $TYPEDFIELDS
"""
struct GeneralizedCirculantBivariateBicycle <: AbstractCSSCode
    """Dimension of cyclic shift matrix `Sₗ` where `x = Sₗ ⊗ Iₘ`"""
    l::Int
    """ Dimension of cyclic shift matrix `Sₘ` where `y = Iₗ ⊗ Sₘ`"""
    m::Int
    """Terms in matrix A, where each tuple is (:x or :y, power)"""
    A::Vector{Tuple{Symbol,Int}}
    """Terms in matrix B, where each tuple is (:x or :y, power)"""
    B::Vector{Tuple{Symbol,Int}}
    
    function GeneralizedCirculantBivariateBicycle(l, m, A, B)
        (l >= 0 && m >= 0) || throw(ArgumentError("l and m must be non-negative"))
        (length(A) >= 1 && length(B) >= 1) || throw(ArgumentError("A and B must each have at least one entry"))
        for (mat, terms) in [(:A, A), (:B, B)]
            for (var, pow) in terms
                var ∈ [:x, :y] || throw(ArgumentError("Matrix $mat contains invalid variable $var (must be :x or :y)"))
                pow >= 0 || throw(ArgumentError("Matrix $mat contains negative power $pow"))
                max_pow = var == :x ? l : m
                pow <= max_pow || throw(ArgumentError("Power $pow in matrix $mat exceeds maximum $max_pow for $var"))
            end
        end
        new(l, m, A, B)
    end
end

function parity_matrix_xz(c::GeneralizedCirculantBivariateBicycle)
    Iₗ = Matrix{Bool}(I, c.l, c.l)
    Iₘ = Matrix{Bool}(I, c.m, c.m)
    xₚ = Dict(i => kron(circshift(Iₗ, (0,i)), Iₘ) for i in 0:c.l)
    yₚ = Dict(i => kron(Iₗ, circshift(Iₘ, (0,i))) for i in 0:c.m)
    A = zeros(Bool, c.l*c.m, c.l*c.m)
    for (var, pow) in c.A
        mat = var == :x ? xₚ[pow] : yₚ[pow]
        A .+= mat
    end
    A = mod.(A, 2)
    B = zeros(Bool, c.l*c.m, c.l*c.m)
    for (var, pow) in c.B
        mat = var == :x ? xₚ[pow] : yₚ[pow]
        B .+= mat
    end
    B = mod.(B, 2)
    Hx = hcat(A, B)
    Hz = hcat(B', A')
    return Hx, Hz
end

code_n(c::GeneralizedCirculantBivariateBicycle) = 2*c.l*c.m

parity_matrix_x(c::GeneralizedCirculantBivariateBicycle) = parity_matrix_xz(c)[1]

parity_matrix_z(c::GeneralizedCirculantBivariateBicycle) = parity_matrix_xz(c)[2]
