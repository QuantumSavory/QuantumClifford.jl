"""
The **La-Cross** code is a quantum LDPC code constructed using the hypergraph product
of two classical seed LDPC codes. It is characterized by its parity check matrix `H`,
which is derived from **circulant** matrices with specific properties. These codes were
introduced in [pecorari2025high](@cite).

The La-Cross code has two families: one for **periodic boundary** conditions and one for
**open boundary** conditions:

```@raw html
<div class="mermaid">
graph TD
    A[La-cross Code Families] --> B[Periodic Boundary]
    A --> C[Open Boundary]

    B -- full_rank = false --> D[⟦2n², 2k², d⟧]
    C -- full_rank = true --> E[⟦❨n - k❩² + n², k², d⟧]
</div>
```

!!! note When `H` is square and circulant (`full_rank=false`), classical checks connect
opposite endpoints of the `length-n` classical code and give rise to a quantum code with
stabilizers connecting opposite array boundaries, i.e. with **periodic boundary conditions**.
On the contrary, **rectangular** parity-check matrices in ``\\mathbb{F}_2^{(n−k)×n}`` give
rise to a quantum code with stabilizers stretching up to the array extent, i.e. with **open
boundary conditions**.

# Cyclic codes and circulant matrices

A **cyclic code** is a linear code in which codewords remain valid under cyclic
shifts. A **circulant matrix** is a square matrix where each row is a cyclic shift
of the first row. When the parity-check matrix `H` is circulant, the code is fully
determined by its first row:

\$\$
H = \\text{circ}(c_0, c_1, \\dots, c_k, 0, \\dots, 0) \\in \\mathbb{F}_2^{n \\times n}.
\$\$

The elements ``c_i`` (for ``i = 0, 1, \\dots, k``) correspond to the coefficients of
a ``degree-k`` polynomial:

\$\$
h(x) = 1 + \\sum_{i=1}^{k} c_i x^i \\in \\mathbb{F}_2[x]/(x^n - 1).
\$\$

This establishes a mapping between ``\\mathbb{F}_2^n`` and the quotient ring
``\\mathbb{F}_2[x]/(x^n - 1)``, where cyclic shifts in ``\\mathbb{F}_2^n`` correspond
to multiplications by `x` in the polynomial ring. Since multiplication by `x` preserves
the ideal structure of ``\\mathbb{F}_2[x]/(x^n - 1)``, cyclic codes correspond to
ideals in this ring. These ideals are in one-to-one correspondence with unitary``mod-2``
divisors of ``x^n - 1`` with a leading coefficient of 1. Consequently, the fundamental
building blocks of cyclic codes correspond to the factorization of ``x^n - 1``.

!!! note For `k = 1`, the generator polynomial `h(x) = 1 + x` defines the **repetition code**.

# Polynomial representation

The first row of a circulant matrix ``H = \\text{circ}(c_0, c_1, c_2, \\dots, c_{n-1})``
can be mapped to the coefficients of a polynomial ``h(x)``. For instance, if the first
row is ``[1, 1, 0, 1]``, the polynomial is: ``h(x) = 1 + x + x^3``. This polynomial-based
representation aids in the analysis and design of cyclic codes.

!!! The **next-to-next-to-nearest neighbor** connectivity implies the use of a *degree-3*
seed polynomial ``h(x) = 1 + x + x^2 + x^3`` in the ring ``\\mathbb{F}_2[x]/(x^n - 1)`` for
a specific code length `n`. Additionally, the condition of low stabilizer weight requires
the polynomial ``1 + x + x^3``.

# [[2n², 2k², d]] La-Cross code

Here is `[[98, 18, 4]]` La-cross code from with `h(x) = 1 + x + x^3`, `n = 7`,
and `k = 3` from Appendix A of [pecorari2025high](@cite).

```jldoctest lacrosseg
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> using QuantumClifford: stab_looks_good

julia> n = 7; k = 3; coeffs = [1,0,1];

julia> c = parity_checks(Lacross(n,k,coeffs,false));

julia> code_n(c), code_k(c)
(98, 18)

julia> stab_looks_good(copy(c), remove_redundant_rows=true)
true
```

# [[(n - k)² + n², k², d]] La-Cross code

Here is `[[65, 9, 4]]` La-cross code from with `h(x) = 1 + x + x^3`, `n = 7`, `k = 3`
and full rank seed *rectangular* circulant matrix from Appendix A of [pecorari2025high](@cite).

```jldoctest lacrosseg
julia> n = 7; k = 3; coeffs = [1,0,1];

julia> c = parity_checks(Lacross(n,k,coeffs,true));

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
    """The degree of the polynomial."""
    k::Int
    """The polynomial representation of the first row of the circulant parity-check matrix."""
    coeffs::Vector{Int}
    """A flag indicating whether to use the full-rank rectangular matrix (true) or the original circulant matrix (false)."""
    full_rank::Bool
    function Lacross(n, k, coeffs, full_rank)
        k <= 0 && throw(ArgumentError("Degree k must be positive."))
        n <= 0 && throw(ArgumentError("Block length must be positive."))
        length(coeffs) != k && throw(ArgumentError("Length of coeffs vector must match degree k."))
        new(n,k,coeffs,full_rank)
    end
end

function iscss(::Type{Lacross})
    return true
end

function parity_checks_xz(c::Lacross)
    first_row = zeros(Int, c.n)
    first_row[1] = 1 # constant term (x⁰)
    for i in 1:c.k
        first_row[i+1] = c.coeffs[i] # coefficients for x¹, x², ..., xᵏ
    end
    H = zeros(Int, c.n, c.n)
    H[1, :] = first_row
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

code_k(c::Lacross) = c.full_rank ? c.k^2 : 2*c.k^2

code_n(c::Lacross) = c.full_rank ? (c.n-c.k)^2+c.n^2 : 2*c.n^2
