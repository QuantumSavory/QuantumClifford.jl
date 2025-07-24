"""
    $TYPEDEF

The **La-cross** code is a quantum LDPC code constructed using the hypergraph product
of two classical seed LDPC codes. It is characterized by its parity check matrix `H`,
which is derived from **circulant** matrices with *specific* properties. These codes were
introduced in [pecorari2025high](@cite).

The La-cross code has two families: one for **periodic boundary** conditions and one for
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

!!! note
    When `H` is square and circulant (`full_rank=false`), classical checks connect opposite
    endpoints of the `length-n` classical code and give rise to a quantum code with stabilizers
    connecting opposite array boundaries, i.e. with **periodic boundary conditions**. On the
    contrary, **rectangular** parity-check matrices in ``\\mathbb{F}_2^{(n−k)×n}`` give rise to
    a quantum code with stabilizers stretching up to the array extent, i.e. with **open boundary
    conditions**.

# Cyclic Codes and Circulant Matrices

A **cyclic code** is a linear code in which codewords remain valid under cyclic
shifts. A **circulant matrix** is a square matrix where each row is a cyclic shift
of the first row. When the parity-check matrix `H` is circulant, the code is fully
determined by its first row:

```math
\\begin{aligned}
H = \\text{circ}(c_0, c_1, \\dots, c_k, 0, \\dots, 0) \\in \\mathbb{F}_2^{n \\times n}.
\\end{aligned}
```

The elements ``c_i`` (for ``i = 0, 1, \\dots, k``) correspond to the coefficients of
a ``degree-k`` polynomial:

```math
\\begin{aligned}
h(x) = 1 + \\sum_{i=1}^{k} c_i x^i \\in \\mathbb{F}_2[x]/(x^n - 1).
\\end{aligned}
```

This establishes a mapping between ``\\mathbb{F}_2^n`` and the quotient ring
``\\mathbb{F}_2[x]/(x^n - 1)``, where cyclic shifts in ``\\mathbb{F}_2^n`` correspond
to multiplications by `x` in the polynomial ring. Since multiplication by `x` preserves
the ideal structure of ``\\mathbb{F}_2[x]/(x^n - 1)``, cyclic codes correspond to
ideals in this ring. These ideals are in one-to-one correspondence with unitary``mod-2``
divisors of ``x^n - 1`` with a leading coefficient of 1. Consequently, the fundamental
building blocks of cyclic codes correspond to the factorization of ``x^n - 1``.

!!! note
    For `k = 1`, the generator polynomial ``h(x) = 1 + x`` defines the **repetition code**.

# Polynomial Representation

The first row of a circulant matrix ``H = \\text{circ}(c_0, c_1, c_2, \\dots, c_{n-1})``
can be mapped to the coefficients of a polynomial ``h(x)``. For instance, if the first
row is ``[1, 1, 0, 1]``, the polynomial is: ``h(x) = 1 + x + x^3``. This polynomial-based
representation aids in the analysis and design of cyclic codes. For our implementation of
La-cross codes, we leverage `Hecke.polynomial_ring` to work directly with polynomial rings
rather than manipulating coefficient arrays explicitly.

!!! The **next-to-next-to-nearest neighbor** connectivity implies the use of a *degree-3*
seed polynomial ``h(x) = 1 + x + x^2 + x^3`` in the ring ``\\mathbb{F}_2[x]/(x^n - 1)`` for
a specific code length `n`. Additionally, the condition of low stabilizer weight requires
the polynomial ``1 + x + x^3``.

## [[2n², 2k², d]] La-cross Code

Here is `[[98, 18, 4]]` La-cross code from with ``h(x) = 1 + x + x^3``, `n = 7`,
and `k = 3` from Appendix A of [pecorari2025high](@cite).

```jldoctest lacrosseg
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> import Hecke: GF, polynomial_ring;

julia> n = 7; k = 3; F = GF(2);

julia> R, x = polynomial_ring(F, "x");

julia> h = 1 + x + x^k;

julia> c = Lacross(n, h, false);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(98, 18, 4)
```

## [[(n - k)² + n², k², d]] La-cross Code

Here is `[[65, 9, 4]]` La-cross code from with ``h(x) = 1 + x + x^3``, `n = 7`, `k = 3`
and full rank seed *rectangular* circulant matrix from Appendix A of [pecorari2025high](@cite).

```jldoctest lacrosseg
julia> c = Lacross(n, h, true);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(65, 9, 4)
```

Here is `[[400, 16, 8]]` La-cross code from with ``h(x) = 1 + x + x^4``, `n = 16`, `k = 4`
from [pecorari2025high](@cite).

```jldoctest lacrosseg
julia> n = 16; k = 4;

julia> R, x = polynomial_ring(F, "x");

julia> h = 1 + x + x^k;

julia> full_rank = true;

julia> c = Lacross(n, h, full_rank);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(400, 16, 8)
```

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/lacross).

### Fields
    $TYPEDFIELDS
"""
struct Lacross <: AbstractCSSCode
    """The block length of the classical seed code."""
    n::Int
    """The seed vector is represented with a degree-`n` polynomial of the form ``h(x) = \\sum_{i=0}^{n-1} c_i x^i``."""
    h::FqPolyRingElem
    """A flag indicating whether to use the full-rank rectangular matrix (`true`) or the original circulant matrix (`false`)."""
    full_rank::Bool
    function Lacross(n, h, full_rank)
        n <= 0 && throw(ArgumentError("Block length must be positive."))
        new(n, h, full_rank)
    end
end

function parity_matrix_xz(c::Lacross)
    F = GF(2)
    R = parent(c.h)
    x = gen(R)
    _, proj = residue_ring(R, R(x)^c.n-1)
    h = proj(c.h)
    lifted_h = lift(h)
    coeffs = Int[lift(ZZ, coeff(lifted_h, i)) for i in 0:c.n-1]
    H = zero_matrix(F, c.n, c.n)
    for i in 1:c.n
        for j in 1:c.n
            H[i, j] = coeffs[mod1(j-i+1, c.n)]
        end
    end
    H = [Int(lift(ZZ, H[i,j])) for i in 1:nrows(H), j in 1:ncols(H)]
    if c.full_rank == true
        k = degree(c.h)
        H = H[1:c.n-k, :]
        hx, hz = hgp(H,H)
        return hx, hz
    else
        hx, hz = hgp(H,H)
        return hx, hz
    end
end

parity_matrix_x(c::Lacross) = parity_matrix_xz(c)[1]

parity_matrix_z(c::Lacross) = parity_matrix_xz(c)[2]

code_k(c::Lacross) = c.full_rank ? degree(c.h)^2 : 2*degree(c.h)^2

code_n(c::Lacross) = c.full_rank ? (c.n-degree(c.h))^2+c.n^2 : 2*c.n^2
