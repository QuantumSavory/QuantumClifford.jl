"""
    $TYPEDEF

The bivariate bicycle (BB) quantum LDPC code was introduced in [bravyi2024high](@cite).
This construction uses identity and cyclic shift matrices because we consider `Iₗ`
as the `l × l` identity matrix and `Sₗ` as the cyclic shift matrix of the same size,
where each row of `Sₗ` has a single '1' at the column `(i + 1) mod l`. We use the
matrices `x = Sₗ ⊗ Iₘ` and `y = Iₗ ⊗ Sₘ`. The BB code is represented by matrices
`A` and `B`, defined as: `A = A₁ + A₂ + A₃` and `B = B₁ + B₂ + B₃`. The addition and
multiplication operations on binary matrices are performed modulo 2. The check matrices
are: `Hx = [A|B]` and `Hz = [B'|A']`. Both `Hx` and `Hz` are `(n/2)×n` matrices.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/qcga).

Here is `[[756, 16, ≤ 34]]` circulant bivariate bicycle code from Table 1 of [bravyi2024high](@cite)
with polynomials ``A = x^3 + y^10 + y^17`` and ``B = y^5 + x^3 + x^19``. 

### Example

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> l, m = 21, 18;

julia> A = [3, 10, 17];

julia> B = [5, 3, 19];

julia> c = CirculantBivariateBicycle(l, m, A, B);

julia> code_n(c), code_k(c)
(756, 16)
```

!!! note
    The Bivariate Bicycle code ``\\mathrm{QC}(A,B)`` is a specific instance of the
    *Lifted Product* construction, where the underlying group is the direct product 
    ``\\mathbb{Z}_\\ell \\times \\mathbb{Z}_m`` (with ``\\mathbb{Z}_j`` denoting the
    cyclic group of order ``j``).

### Fields
    $TYPEDFIELDS
"""
struct CirculantBivariateBicycle <: AbstractCSSCode
    """Dimension of cyclic shift matrix ``S_\\ell`` where ``x = S_\\ell \\otimes I_m``."""
    l::Int
    """Dimension of cyclic shift matrix ``S_m`` where ``y = I_\\ell \\otimes S_m``."""
    m::Int
    """Specifies ``[A_1, A_2, A_3]`` as powers of ``x`` or ``y\``."""
    A::Vector{Int}
    """Specifies ``[B_1, B_2, B_3]`` as powers of ``x`` or ``y``."""
    B::Vector{Int}
    function CirculantBivariateBicycle(l,m,A,B)
        (l >= 0 && m >= 0) || throw(ArgumentError("l and m must be non-negative"))
        (length(A) == 3 && length(B) == 3) || throw(ArgumentError("A and B must each have exactly 3 entries"))
        (all(x -> x >= 0, A) && all(x -> x >= 0, B)) || throw(ArgumentError("A and B must contain only non-negative integers"))
        (all(x -> x in 0:max(l, m), A) && all(x -> x in 0:max(l, m), B)) || throw(ArgumentError("Each element in A and B must be in the range [0, $(max(l, m))]."))
        new(l,m,A,B)
    end
end

function parity_matrix_xz(c::CirculantBivariateBicycle)
    a₁,a₂,a₃ = c.A[1],c.A[2],c.A[3]
    b₁,b₂,b₃ = c.B[1],c.B[2],c.B[3]
    Iₗ  = Matrix{Bool}(LinearAlgebra.I,c.l,c.l)
    Iₘ = Matrix{Bool}(LinearAlgebra.I,c.m,c.m)
    x  = Dict{Bool, Matrix{Bool}}()
    y  = Dict{Bool, Matrix{Bool}}()
    x  = Dict(i => kron(circshift(Iₗ,(0,i)),Iₘ) for i in 0:(c.l))
    y  = Dict(i => kron(Iₗ,circshift(Iₘ,(0,i))) for i in 0:(c.m))
    A  = mod.(x[a₁]+y[a₂]+y[a₃],2)
    B  = mod.(y[b₁]+x[b₂]+x[b₃],2)
    Hx = hcat(A,B)
    Hz = hcat(B',A')
    return Hx, Hz
end

code_n(c::CirculantBivariateBicycle) = 2*c.l*c.m

parity_matrix_x(c::CirculantBivariateBicycle) = parity_matrix_xz(c)[1]

parity_matrix_z(c::CirculantBivariateBicycle) = parity_matrix_xz(c)[2]
