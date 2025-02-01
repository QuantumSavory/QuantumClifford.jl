"""
A La-cross code is quantum LDPC code constructed using the hypergraph product of two
classical LDPC codes. The LaCross LDPC code is characterized by its parity-check matrix,
which is derived from circulant matrices with specific properties.

# Properties of Circulant Matrices

A circulant matrix is one where each row is a cyclic shift of the first row. The
parity-check matrix H can be described as:

\$\$
H = \\text{circ}(c_0, c_1, c_2, \\dots, c_{n-1}) \\in \\mathbb{F}_2^{n \\times n}.
\$\$

Each element ``c_i`` (for ``i = 0, 1, \\dots, n-1``) corresponds to a polynomial of
degree n-1:

\$\$
h(x) = c_0 + c_1 x + c_2 x^2 + \\dots + c_{n-1} x^{n-1}.
\$\$


This polynomial belongs to the ring ``\\mathbb{F}_2[x]/(x^n - 1)``, which consists of
polynomials modulo ``x^n - 1``. In this formulation, cyclic shifts in ``\\mathbb{F}_2^n``
correspond to multiplications by ``x`` in ``\\mathbb{F}_2[x]/(x^n - 1)``, making this
representation particularly useful in coding theory.

# Polynomial Representation

The first row of a circulant matrix ``H = \\text{circ}(c_0, c_1, c_2, \\dots, c_{n-1})``
can be mapped to the coefficients of a polynomial ``h(x)``. For instance, if the first
row is ``[1, 1, 0, 1]``, the polynomial is: ``h(x) = 1 + x + x^3``. This polynomial-based
representation aids in the analysis and design of cyclic codes.

# Example

An `[[98, 18, 4]]` La-cross code from with `h(x) = h(x) = 1 + x + x^3` and `n = 7`
from [pecorari2025high](@cite).

```jldoctest lacrosseg
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> using QuantumClifford: stab_looks_good

julia> n = 7; polynomial = [1,1,0,1];

julia> c = parity_checks(Lacross(n,polynomial,false));

julia> code_n(c), code_k(c)
(98, 18)

julia> stab_looks_good(copy(c), remove_redundant_rows=true)
true
```

An `[[65, 9, 4]]` La-cross code from with `h(x) = h(x) = 1 + x + x^3`, `n = 7` and
full rank seed circulant matrix from [pecorari2025high](@cite).

```jldoctest lacrosseg
julia> n = 7; polynomial = [1,1,0,1];

julia> c = parity_checks(Lacross(n,polynomial,true));

julia> code_n(c), code_k(c)
(65, 9)

julia> stab_looks_good(copy(c), remove_redundant_rows=true)
true
```

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/lacross)
"""
struct Lacross <: AbstractECC
    """The block length of the classical seed code"""
    n::Int
    """The polynomial representation of the first row of the circulant parity-check matrix."""
    pattern::Vector{Int}
    """A flag indicating whether to use the full-rank rectangular matrix (true) or the original circulant matrix (false)."""
    full_rank::Bool
end

function iscss(::Type{Lacross})
    return true
end

function parity_checks_xz(c::Lacross)
    first_r = zeros(Int, c.n)
    first_r[1:length(c.pattern)] = c.pattern
    H = zeros(Int, c.n, c.n)
    H[1, :] = first_r
    for i in 2:c.n
        H[i, :] = circshift(H[i-1, :], 1)
    end
    if c.full_rank == true
        rank = QuantumClifford.gf2_row_echelon_with_pivots!(H)[2]
        H = H[1:rank, :]
        hx, hz = hgp(H,H)
        return hx, hz
    else
        hx, hz = hgp(H,H)
        return hx, hz
    end
end

parity_checks_x(c::Lacross) = parity_checks_xz(c)[1]

parity_checks_z(c::Lacross) = parity_checks_xz(c)[2]

parity_checks(c::Lacross) = parity_checks(CSS(parity_checks_xz(c)...))
