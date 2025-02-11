using Random: AbstractRNG, GLOBAL_RNG, randperm

"""
The **(n, m, r)-Structured quantum LDPC codes** code is constructed using
the hypergraph product of two classical seed **structured LDPC** codes.

The classical structured LDPC were introduced in [tillich2006minimum](@cite).
Arnault et. al. showed ([arnault2025upperboundsminimumdistance](@cite)) that
the minimum distance of structured binary LDPC codes with parity-check matrices
of the form [C|M], where C is circulant and M has fixed column weight ``r \\geq 3``,
is in ``O(n^{\\frac{r-2}{r-1} + \\epsilon})``, improving the previous bound of
``O(n^{\\frac{r-1}{r}}``.

The classical **structured LDPC** codes are defined by a specific structure in their
parity-check matrix `H`, which is composed of a circulant matrix `C` and a binary
matrix `M`. The ensemble of these codes is referred to as *(n, m, r)-structured LDPC*
codes, where `n` is the block length of the code, `m` is the number of parity-check
equations, `r` is the weight of the columns in the matrix `M`.

## Parity-Check Matrix Structure

The parity-check matrix `H` is an `m × n` binary matrix of the form: `H = [C | M]`
where: `C`  is an  ``m \\times m``  **circulant matrix**. `M` is an  ``m \\times (n - m)``
binary matrix with specific properties.

### Circulant Matrix C

The circulant matrix `C` is defined as:

\$\$
C = \\begin{pmatrix}
1 &        &        &        & 1 \\\\
1 & 1      &        &        &   \\\\
  & 1      & \\ddots &        &   \\\\
  &        & \\ddots & 1      &   \\\\
  &        &        & 1      & 1
\\end{pmatrix}
\$\$

This matrix has the property that each row is a cyclic shift of the previous
row. The first row contains two consecutive 1s, and the rest of the entries are 0s.

### Binary Matrix M

The matrix `M` satisfies the following conditions:
1. **No Zero Rows**: Every row of `M` contains at least one non-zero element.
2. **Column Weight**: Every column of `M` has a constant weight `r`, where ``r \\geq 3``.

## Existence Conditions

The family of classical *(n, m, r)-structured LDPC* codes exists if and only if the
following conditions are satisfied:

\$\$
m \\geq r \\quad \\text{and} \\quad (n - m)r \\geq m
\$\$.
These conditions ensure that the matrix `M` can be constructed with the required properties.
"""
struct StructuredQLDPC <: AbstractECC
    """The block length of the classical seed code"""
    n::Int
    """The number of check nodes (rows in the parity-check matrix H)"""
    m::Int
    """The column weight parameter for matrix M (each column of M must have exactly r ones)"""
    r::Int
    function StructuredQLDPC(n, m, r)
        m < r || (n - m) * r < m && throw(ArgumentError(("Conditions for the existence of M are not satisfied.")))
        new(n, m, r)
    end
end

function iscss(::Type{ StructuredQLDPC})
    return true
end

"""Function to create a circulant matrix C of size m x m"""
function _create_circulant_matrix(m::Int)
    first_row = zeros(Int, m)
    first_row[1] = 1
    first_row[end] = 1
    C = zeros(Int, m, m)
    for i in 1:m
        C[i, :] = circshift(first_row, i-1)
    end
    return C
end

"""create the matrix M of size m x (n - m) with column weight r"""
function _create_matrix_M_deterministic(m::Int, n::Int, r::Int)
    M = zeros(Int, m, n - m)
    # Fill M such that each column has exactly r ones
    for col in 1:(n - m)
        # Deterministically select r distinct rows to place ones
        rows = mod1.((col-1)*r+1:col*r, m)
        M[rows, col] .= 1
    end
    # Ensure no row is all zeros
    for row in 1:m
        if all(M[row, :] .== 0)
            # Deterministically select the first column to set this row to 1
            M[row, 1] = 1
        end
    end
    return M
end

function parity_checks_xz(c:: StructuredQLDPC)
    C = _create_circulant_matrix(c.m)
    M = _create_matrix_M_deterministic(c.m, c.n, c.r)
    # The parity-check matrix H = [C | M]
    H = hcat(C, M)
    hx, hz = hgp(H,H)
    return hx, hz
end

parity_checks_x(c:: StructuredQLDPC) = parity_checks_xz(c)[1]

parity_checks_z(c:: StructuredQLDPC) = parity_checks_xz(c)[2]

parity_checks(c:: StructuredQLDPC) = parity_checks(CSS(parity_checks_xz(c)...))

"""
The **random (n, m, r)-Structured quantum LDPC codes** code is constructed
using the hypergraph product of two classical seed **random structured LDPC** codes.
"""
function random_structured_qldpc_code end

function _create_matrix_M_random(rng::AbstractRNG, m::Int, n::Int, r::Int)
    M = zeros(Int, m, n - m)
    for col in 1:(n - m)
        # Randomly select r distinct rows to place ones
        rows = randperm(rng, m)[1:r]
        M[rows, col] .= 1
    end
    # Ensure no row is all zeros
    for row in 1:m
        if all(M[row, :] .== 0)
            # Randomly select a column and set this row to 1
            col = rand(rng, 1:(n - m))
            M[row, col] = 1
        end
    end
    return M
end

function _construct_parity_check_matrix(rng::AbstractRNG, n::Int, m::Int, r::Int)
    C = _create_circulant_matrix(m)
    M = _create_matrix_M_random(rng, m, n, r)
    # The parity-check matrix H = [C | M]
    H = hcat(C, M)
    return H
end

function random_structured_qldpc_code(rng::AbstractRNG, n::Int, m::Int, r::Int)
    H = _construct_parity_check_matrix(rng, n, m, r)
    hx, hz = hgp(H, H)
    Stabilizer(CSS(hx, hz))
end

function random_structured_qldpc_code(n::Int, m::Int, r::Int)
    H = _construct_parity_check_matrix(GLOBAL_RNG, n, m, r)
    hx, hz = hgp(H, H)
    Stabilizer(CSS(hx, hz))
end
