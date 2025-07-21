"""
    $TYPEDEF

The `[[n² + m²,(n - rank([C ∣ M]))² + (m − rank([C ∣ M]ᵀ))², d]]` **quantum
Tillich Zémor code** is a novel quantum LDPC code is constructed using the
hypergraph product of two classical seed **(n, m, r)-Structured LDPC** codes.

# Structured LDPC

The classical structured LDPC were introduced in [tillich2006minimum](@cite).
Arnault et. al. showed ([arnault2025upperboundsminimumdistance](@cite)) that
the minimum distance of structured binary LDPC codes with parity-check matrices
of the form ``[C \\mid M]``, where `C` is circulant and `M` has fixed column
weight ``r \\geq 3``, is in ``O(n^{\\frac{r-2}{r-1} + \\epsilon})``, improving
the previous bound of ``O(n^{\\frac{r-1}{r}}``.

The classical **structured LDPC** codes are defined by a specific structure in their
parity-check matrix `H`, which is composed of a circulant matrix `C` and a binary
matrix `M`. The ensemble of these codes is referred to as *(n, m, r)-structured LDPC*
codes, where `n` is the block length of the code, `m` is the number of parity-check
equations, `r` is the weight of the columns in the matrix `M`.

## Parity-Check Matrix Structure

The parity-check matrix `H` is an `m × n` binary matrix of the form: ``[C \\mid M]``
where: `C`  is an  ``m \\times m``  **circulant matrix**. `M` is an  ``m \\times (n - m)``
binary matrix with specific properties.

### Circulant Matrix C

The circulant matrix `C` is defined as:

```math
\\begin{aligned}
C = \\begin{pmatrix}
1 &        &        &        & 1 \\\\
1 & 1      &        &        &   \\\\
  & 1      & \\ddots &        &   \\\\
  &        & \\ddots & 1      &   \\\\
  &        &        & 1      & 1
\\end{pmatrix}
\\end{aligned}
```

This matrix has the property that each row is a cyclic shift of the previous
row. The first row contains two consecutive 1s, and the rest of the entries are 0s.

### Binary Matrix M

The matrix `M` satisfies the following conditions:
- **No Zero Rows**: Every row of `M` contains at least one non-zero element.
- **Column Weight**: Every column of `M` has a constant weight `r`, where ``r \\geq 3``.

## Existence Conditions

The family of classical *(n, m, r)-structured LDPC* codes exists if and only if the
following conditions are satisfied:

```math
\\begin{aligned}
m \\geq r \\quad \\text{and} \\quad (n - m)r \\geq m
\\end{aligned}
```

These conditions ensure that the matrix `M` can be constructed with the required properties.


# Quantum Structured LDPC codes

We introduce a novel construction of quantum LDPC codes, inspired by classical
structured LDPC seed codes. Leveraging the advantages of classical structured LDPC
codes—namely, their linear-time encoding and sub-linear scaling of minimum distance
—we present the `[[n² + m²,(n - rank([C ∣ M]))² + (m − rank([C ∣ M]ᵀ))², d]]` **quantum
Tillich Zémor** code. Our approach employs structure parity-check matrices of the form
``H = [C \\mid M]``, where `C` serves as the circulant core matrix and `M` is meticulously
designed to maintain crucial code properties while enabling efficient encoding. This QECC
emphasizes the potential of using structured classical LDPC codes as a robust framework for
develop ing scalable and effective QECCs, particularly the quantum Tillich-Zémor code.

# Examples

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> n, m, r = 4, 3, 3;

julia> c = parity_checks(TillichZemor(n, m, r))
+ X_______X___X___XX_______
+ _X_______X___X___XX______
+ __X_______X___X_X_X______
+ ___X_______X___XXXX______
+ X___X_______X______XX____
+ _X___X_______X______XX___
+ __X___X_______X____X_X___
+ ___X___X_______X___XXX___
+ ____X___X___X_________XX_
+ _____X___X___X_________XX
+ ______X___X___X_______X_X
+ _______X___X___X______XXX
+ Z_ZZ____________Z__Z_____
+ ZZ_Z_____________Z__Z____
+ _ZZZ______________Z__Z___
+ ____Z_ZZ___________Z__Z__
+ ____ZZ_Z____________Z__Z_
+ _____ZZZ_____________Z__Z
+ ________Z_ZZ____Z_____Z__
+ ________ZZ_Z_____Z_____Z_
+ _________ZZZ______Z_____Z
+ ____________Z_ZZZ__Z__Z__
+ ____________ZZ_Z_Z__Z__Z_
+ _____________ZZZ__Z__Z__Z

julia> code_n(c), code_k(c)
(25, 1)

julia> n, m, r = 100, 40, 40;

julia> c = TillichZemor(n, m, r);

julia> code_n(c), code_k(c)
(11600, 3722)
```

### Fields
    $TYPEDFIELDS
"""
struct TillichZemor{M} <: AbstractCSSCode where {M <: Union{Nothing,Tuple{AbstractMatrix,AbstractMatrix}}}
    """The block length of the classical seed code."""
    n::Int
    """The number of check nodes (rows in the parity-check matrix `H`)."""
    m::Int
    """The column weight parameter for matrix `M` (each column of M must have exactly `r` ones)."""
    r::Int
    """The `X`-type and `Z`-type parity check matrices generated via the hypergraph product of the classical
    parity check matrix `H = [C | M]`. For randomized constructions via `random_TillichZemor_code`, these
    store the matrices from the randomly generated seed code."""
    matrices::M

    function TillichZemor(n::Int, m::Int, r::Int, matrices::M=nothing) where {M <: Union{Nothing,Tuple{AbstractMatrix,AbstractMatrix}}}
        (m ≥ r && (n - m)*r ≥ m) || throw(ArgumentError(("Conditions for the existence of `M` in `H = [C | M]` are not satisfied.")))
        new{M}(n, m, r, matrices)
    end
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

"""create the matrix M of size m x (n - m) with column weight r."""
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

function parity_matrix_xz(c::TillichZemor{Nothing})::Tuple{AbstractMatrix,AbstractMatrix}
    C = _create_circulant_matrix(c.m)
    M = _create_matrix_M_deterministic(c.m, c.n, c.r)
    # The parity-check matrix H = [C | M]
    H = hcat(C, M)
    hx, hz = hgp(H,H)
    return hx, hz
end

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
    (m ≥ r && (n - m)*r ≥ m) || throw(ArgumentError(("Conditions for the existence of `M` in `H = [C | M]` are not satisfied.")))
    C = _create_circulant_matrix(m)
    M = _create_matrix_M_random(rng, m, n, r)
    # The parity-check matrix H = [C | M]
    H = hcat(C, M)
    return H
end

"""
The **random Tillich Zémor code** is a quantum LDPC code constructed using the
hypergraph product of two classical seed **(n, m, r)-Structured LDPC** codes.
"""
function random_TillichZemor_code(rng::AbstractRNG, n::Int, m::Int, r::Int)
    H = _construct_parity_check_matrix(rng, n, m, r)
    hx, hz = hgp(H, H)
    TillichZemor(n, m, r, (hx, hz))
end

function random_TillichZemor_code(n::Int, m::Int, r::Int)
    H = _construct_parity_check_matrix(GLOBAL_RNG, n, m, r)
    hx, hz = hgp(H, H)
    TillichZemor(n, m, r, (hx, hz))
end

function parity_matrix_xz(c::TillichZemor{<:Tuple})::Tuple{AbstractMatrix,AbstractMatrix}
    return c.matrices
end

parity_matrix_x(c::TillichZemor) = parity_matrix_xz(c)[1]

parity_matrix_z(c::TillichZemor) = parity_matrix_xz(c)[2]

code_n(c::TillichZemor) = c.n^2 + c.m^2
