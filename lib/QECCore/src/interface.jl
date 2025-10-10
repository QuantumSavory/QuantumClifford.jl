"""
    AbstractECC

Abstract type for error correction code.
"""
abstract type AbstractECC end

"""
    AbstractQECC <: AbstractECC

Abstract type for quantum error correction code.
"""
abstract type AbstractQECC <: AbstractECC end

"""
    AbstractCECC <: AbstractECC

Abstract type for classical error correction code.
"""
abstract type AbstractCECC <: AbstractECC end

"""
    parity_matrix(c::AbstractECC)

The parity check matrix of a error correction code in the form of `(X|Z)`. The size of the matrix is `(code_s, 2*code_n)`. `code_s` is the number of stabilizers, and `code_n` is the number of physical qubits. Each row of the matrix is a stabilizer. The first `code_n` columns represent whether this stabilizer contains a X operator on the physical qubit, and the last `code_n` columns represent whether this stabilizer contains a Z operator on the physical qubit.

See also: [`parity_matrix_x`](@ref) and [`parity_matrix_z`](@ref)
"""
function parity_matrix end

"""
    AbstractCSSCode <: AbstractQECC

Abstract type for Calderbank-Shor-Steane (CSS) code.
"""
abstract type AbstractCSSCode <: AbstractQECC end

"""
    parity_matrix_x(c::AbstractCSSCode)

Parity check boolean matrix of a code (only the X entries in the tableau, i.e. the checks for Z errors).
Only CSS codes have this method.

See also: [`parity_matrix`](@ref) and [`parity_matrix_z`](@ref)
"""
function parity_matrix_x end

"""
    parity_matrix_z(c::AbstractCSSCode)

Parity check boolean matrix of a code (only the Z entries in the tableau, i.e. the checks for X errors).
Only CSS codes have this method.

See also: [`parity_matrix`](@ref) and [`parity_matrix_x`](@ref)
"""
function parity_matrix_z end


"""
    code_n(c::AbstractECC)

The number of physical qubits in a error correction code.

See also: [`code_k`](@ref) and [`code_s`](@ref)
"""
code_n(c::AbstractQECC) = nqubits(parity_matrix(c))
code_n(c::AbstractCECC) = nbits(parity_matrix(c))
nqubits(pm::AbstractMatrix{Bool}) = size(pm, 2) ÷ 2
nbits(pm::AbstractMatrix{Bool}) = size(pm, 2)

"""
    code_s(c::AbstractECC)

The number of stabilizers in a error correction code. They might not be all linearly independent, thus `code_s >= code_n-code_k`. For the number of linearly independent checks you can use `LinearAlgebra.rank`.

See also: [`code_n`](@ref) and [`code_k`](@ref)
"""
function code_s end
code_s(c::AbstractQECC) = nstabilizers(parity_matrix(c))
nstabilizers(pm::AbstractMatrix{Bool}) = size(pm, 1)

"""
    code_k(c::AbstractECC)

The number of logical qubits in a error correction code.

See also: [`code_n`](@ref) and [`code_s`](@ref)
"""
function code_k end

"""
    rate(c::AbstractECC)

The rate of a error correction code.

See also: [`code_n`](@ref) and [`code_k`](@ref)
"""
function rate(c)
    rate = code_k(c)//code_n(c)
    return rate
end

"""
    distance(c::AbstractECC)

The code distance of a error correction code.

See also: [`code_n`](@ref) and [`code_k`](@ref)
"""
function distance end

"""
    AbstractDistanceAlg

Abstract type representing algorithms for computing
the minimum distance of quantum error correction codes.
"""
abstract type AbstractDistanceAlg end

"""
    metacheck_matrix_x(c::AbstractCSSCode)

Returns the `X`-metacheck matrix (``M_X = \\partial_{i-1}``) for a CSS code defined
by a chain complex of length `l ≥ 4`, where qubits are placed on `i`-cells with
``1 < i < l−1``. This matrix acts on `X`-syndromes (measurement outcomes of `Z`-type
stabilizers) and enforces the constraint ``M_Xs_X = 0``, ensuring syndromes are valid
codewords of a classical *metacode*.

!!! note
    For an introduction to chain complexes in quantum error correction and the role
    of metachecks in single-shot QEC, see the documentation for [`metacheck_matrix`](@ref).

### Example: 4D surface code

In the 5-term chain complex used for the 4D surface code:

```math
\\begin{aligned}
C_4 \\xrightarrow{\\partial_4} C_3 \\xrightarrow{\\partial_3} C_2 \\xrightarrow{\\partial_2} C_1 \\xrightarrow{\\partial_1} C_0
\\end{aligned}
```

the metacheck matrix ``M_X = \\partial_1`` satisfies the following:

- Acts on `X`-syndromes: ``s_X \\in C_1``.
- It enforces ``M_Xs_X = 0``, i.e. only valid syndromes lie in ``\\ker M_X``.
- It satisfies the boundary condition ``M_XH_X = 0`` (i.e., ``\\partial_1 \\partial_2 = 0``).

Only CSS codes built using chain complexes and homology have this method.

See also: [`metacheck_matrix_z`](@ref), [`metacheck_matrix`](@ref), [`parity_matrix_x`](@ref)
"""
function metacheck_matrix_x end

"""
    metacheck_matrix_z(c::AbstractCSSCode)

Returns the `Z`-metacheck matrix (``M_Z = \\partial_{i+2}^\\top``) for a CSS code
defined by a chain complex of length `l ≥ 4`, where qubits are placed on `i`-cells
(`1 < i < l−1`).  This matrix validates `Z`-syndromes (measurement outcomes of `X`-type
stabilizers) by enforcing ``M_Zs_Z = 0``, ensuring syndromes are codewords of a classical *metacode*.

!!! note
    For an introduction to chain complexes in quantum error correction and the role
    of metachecks in single-shot QEC, see the documentation for [`metacheck_matrix`](@ref).

### Example: 4D Surface Code

In the `5`-term chain complex used for the `4D` surface code:

```math
\\begin{aligned}
C_4 \\xrightarrow{\\partial_4} C_3 \\xrightarrow{\\partial_3} C_2 \\xrightarrow{\\partial_2} C_1 \\xrightarrow{\\partial_1} C_0
\\end{aligned}
```

the metacheck matrix ``M_Z = \\partial_4`` satisfies the following:

- Acts on Z-syndromes: ``s_Z \\in C_3``,
- It enforces ``M_Z\\cdot_Z = 0``, i.e. only valid syndromes lie in ``\\ker M_Z``.
- It satisfies the boundary condition ``H_Z^\\top \\cdot M_Z^\\top = 0`` (i.e., \\partial_3 \\partial_4 = 0``,

Only CSS codes built using chain complexes and homology have this method.

See also: [`metacheck_matrix_x`](@ref), [`metacheck_matrix`](@ref), [`parity_matrix_z`](@ref)
"""
function metacheck_matrix_z end

"""
    metacheck_matrix(c::AbstractCSSCode)

Return the `X`- and `Z`-metacheck matrices for CSS codes enabling **single-shot
quantum error correction** — a fault-tolerant scheme that corrects both data and
measurement errors using **one** round of syndrome measurements.

### Single-Shot QEC

Single-shot QEC enables both physical errors on qubits and errors in syndrome measurements
to be detected and corrected using only a *single round of noisy measurements*, without
requiring repeated measurement rounds. A layer of redundancy is added to the measurement
process itself. This redundancy is captured by *metasyndromes*: linear constraints that the
noisy syndrome outcomes must satisfy if no measurement error has occurred. When violated, they
indicate faults in the syndrome extraction layer itself.

Traditional QEC combats measurement faults by repeating stabilizer measurements. In contrast,
**single-shot QEC** uses *spatial redundancy* via **metachecks** — extra linear constraints
on syndrome outcomes ("checks of checks") to detect and correct measurement errors immediately.

### Confinement and Single-Shot Decoding

Metachecks enable single-shot decoding by providing a metacode for syndrome repair,
but their role is best understood through the broader property of **confinement**:

- A code has ``*(t,f)*-confinement if, for all errors *E* with ``|E|_{\\text{red}} ``\\leq t``,
the syndrome weight satisfies ``f(|\\sigma(E)|) \\geq |E|_{\\text{red}}``, where
``f: \\mathbb{Z} \\to \\mathbb{R}`` is an increasing function. This bounds the physical
error weight by a function of the syndrome weight.

Codes with metachecks (e.g., `D`-dimensional surface and toric codes) exhibit confinement
because ``M_X/M_Z`` constrain syndromes to a metacode, but confinement can exist *without*
metachecks (e.g., quantum expander codes). Thus, while metachecks are sufficient for single-shot
QEC (via syndrome repair), they are not strictly necessary.

### Repair-Syndrome Decoding

To correct errors in CSS codes, a **two-stage** decoder can be employed when given a noisy
syndrome measurement `z′`. This method separately addresses data qubit errors (e.g., `Z`-errors)
and syndrome measurement errors (e.g., faulty `X`-stabilizer measurements). The same approach
applies symmetrically for `X`-errors and `Z`-stabilizer measurements.

#### Stage I: Syndrome Repair via Metachecks

- The metacheck matrix ``M_X`` computes the metasyndrome ``s = M_Xz'``, identifying inconsistencies caused by measurement errors.
- A classical decoder ``f^1_d : \\mathbb{F}^{n_{i-2}}_2 \\rightarrow \\mathbb{F}^{n_{i-1}}_2`` estimates the noiseless syndrome ``z`` from ``s``, effectively "repairing" the syndrome.

#### Stage II: Data Qubit Correction

A second decoder ``f^2_d : \\mathbb{F}^{n_{i-1}}_2 \\rightarrow \\mathbb{F}^n_2``​ uses the
corrected syndrome `z` to compute a noise vector `n` such that ``H_{X}n = z``, determining
the most likely data qubit errors.

The correction can fail in two ways:
- **Invalid Syndrome**: The repaired `z` lies outside the valid syndrome space ``im(H_X​``).
- **Logical Error**: The correction `n` corresponds to a nontrivial logical operator (i.e., ``n \\in \\ker(H_X)``).

To mitigate the first failure (i.e. invalid syndrome), we can modify the metacheck matrix to

```math
\\begin{aligned}
M' = \\begin{pmatrix}
M_X \\\\
L_M
\\end{pmatrix}
\\end{aligned}
```

where ``L_M`` spans the cohomology group ``H^{i-1}``. However, if ``L_M`` is non-sparse
(common in topological codes), decoding efficiency may suffer.

!!! note
    The modified metacheck matrix ``M'`` is only employed when the initial decoding
    attempt yields an invalid syndrome (``z \\notin \\text{im}(H_X)``). This approach
    helps maintain decoder efficiency, particularly for topological codes where ``L_M`` is
    typically non-sparse. However, the overall performance of two-stage decoders remains
    fundamentally constrained by the metacode's threshold, often resulting in suboptimal
    error correction capability.

### Single-Stage Decoding

Single-stage decoding provides a unified framework for correcting both data qubit errors
(`e`) and syndrome measurement errors (`s_e`) simultaneously. Given an observed syndrome
``s = H_X e + s_e``, where `H_X` is the `X`-stabilizer matrix, the decoder seeks the most
probable error configuration ``e' = \\begin{pmatrix}e \\\\ s_e\\end{pmatrix}`` that satisfies
the extended parity-check equation:

```math
\\begin{aligned}
H' e' = s \\quad \\text{where} \\quad H' = \\begin{pmatrix}H_X & I_r\\end{pmatrix}
\\end{aligned}
```

The Tanner graph ``T(H')`` for this system is constructed by augmenting the original Tanner
graph ``T(H_X)`` with additional variable nodes ``\\{v_i^m\\}_{i=1}^r`` representing potential
measurement errors. Each check node ``f_i`` gains a corresponding edge ``(v_i^m, f_i)``, creating
a structure where measurement errors appear explicitly in the decoding graph.

The decoding framework can be further enhanced by explicitly incorporating *metachecks* through
the extended matrix:

```math
\\begin{aligned}
\\begin{equation}
H_M = \begin{pmatrix}
H_X & I_r \\\\
0 & M
\\end{pmatrix}
\\end{equation}
\\end{aligned}
```

where ``M`` is the metacheck matrix. Though these metachecks are implicitly present as linear
combinations in ``T(H')``, their explicit inclusion significantly improves decoder performance.

The single-stage decoding approach offers several key advantages over two-stage methods:

- elimination of metacode failures (since ``s + s_e \\in \\text{im}(H_X)`` by construction)
- avoidance of non-sparse ``L_M`` matrices that degrade decoder performance
- improved thresholds for topological codes.

### Chain Complexes and ``\\mathbb{F_2}`` Homology

A chain complex of length `l` is a sequence of vector spaces connected by boundary maps:

```math
\\begin{aligned}
\\{0\\} \\xrightarrow{\\partial_{l+1}} C_l \\xrightarrow{\\partial_l} C_{l-1} \\xrightarrow{\\partial_{l-1}} \\cdots \\xrightarrow{\\partial_1} C_0 \\xrightarrow{\\partial_0} \\{0\\}
\\end{aligned}
```

where

- Each ``C_i`` is called an *i-cell*.
- The image of ``\\partial_{i+1}``, denoted ``\\mathrm{im}\\partial_{i+1}``, consists of *i-boundaries*.
- The kernel of ``\\partial_i``, denoted ``\\ker\\partial_i``, consists of *i-cycles*.

The boundary maps satisfy the constraint:

```math
\\begin{aligned}
\\partial_i \\circ \\partial_{i+1} = 0 \\quad \\text{for all } i \\in \\{0, \\dots, l\\}
\\end{aligned}
```

Because ``\\partial_i \\circ \\partial_{i+1} = 0``, every boundary is a cycle:

```math
\\begin{aligned}
\\mathrm{im}\\partial_{i+1} \\subseteq \\ker\\partial_i
\\end{aligned}
```

The **i-th homology group** measures the difference between cycles and boundaries:

```math
\\begin{aligned}
H_i = \\frac{\\ker\\partial_i}{\\mathrm{im}\\partial_{i+1}}
\\end{aligned}
```

Associated with a chain complex is a **cochain complex** with *coboundary operators*
``\\delta^i: C^i \\to C^{i+1}``, typically defined as the transpose (or dual) of the boundary maps:

```math
\\begin{aligned}
\\{0\\} \\xrightarrow{\\delta^{-1}} C^0 \\xrightarrow{\\delta^0} C^1 \\xrightarrow{\\delta^1} \\cdots \\xrightarrow{\\delta^{l-1}} C^l \\xrightarrow{\\delta^l} \\{0\\}
\\end{aligned}
```

where

- ``\\ker\\delta^i`` consists of *i-cocycles*.
- ``\\mathrm{im}\\delta^{i-1}`` consists of *i-coboundaries*.

The **i-th cohomology group** is:

```math
\\begin{aligned}
H^i = \\frac{\\ker\\delta^i}{\\mathrm{im}\\delta^{i-1}}
\\end{aligned}
```

### CSS codes using Homological Algebra

Quantum CSS codes can be described using the framework of [chain complexes](https://en.wikipedia.org/wiki/Chain_complex).

For a chain complex of length `l ≥ 4` , where qubits are placed on `i`-cells (`C_i`) with (`1 < i < l−1`):

```math
\\begin{aligned}
C_{l-1} \\xrightarrow{\\partial_{l-1}} \\cdots \\xrightarrow{\\partial_{i+2}} C_{i+1} \\xrightarrow{\\partial_{i+1}} C_i \\xrightarrow{\\partial_i} C_{i-1} \\xrightarrow{\\partial_{i-1}} \\cdots \\xrightarrow{\\partial_1} C_0
\\end{aligned}
```

where

- **X-stabilizers** are given by the boundary map ``H_X = \\partial_i: C_i → C_{i-1}``.
- **Z-stabilizers** are given by the coboundary map ``H_Z = \\partial_{i+1}^T: C_i → C_{i+1}``.
- **X-metachecks** are defined as ``M_X = \\partial_{i-1}: C_{i-1} → C_{i-2}``.
- **Z-metachecks** are defined as ``M_Z = \\partial_{i+2}^T: C_{i+2} → C_{i+1}``.

The boundary conditions ``\\partial_{i-1} \\partial_i = 0`` (i.e., ``M_X H_X = 0``) guarantee that valid syndromes
(`im H_X`) lie in `ker M_X`.

Invalid syndromes in `ker M_X \\setminus im H_X` belong to the `(i−1)`-th homology group
``H_{i-1} = \\ker \\partial_{i-1} / \\mathrm{im} \\partial_i``, while invalid `Z`-syndromes in
`ker M_Z \\setminus im H_Z` belong to the `(i+1)`-th cohomology group.

!!! note A code can be designed to incorporate syndromes within a metacode by employing
    a chain complex of minimum length three—sufficient for encoding either `X` or `Z` syndromes.
    If the goal is to include both `X` and `Z` syndromes in the metacode, the chain complex must
    extend to at least length four.

### Metachecks in Higher-Dimensional Complexes

In D-dimensional codes, such as the `4D` surface code, we consider a `5`-term chain complex:

```math
\\begin{aligned}
C_4 \\xrightarrow{\\partial_4} C_3 \\xrightarrow{\\partial_3} C_2 \\xrightarrow{\\partial_2} C_1 \\xrightarrow{\\partial_1} C_0
\\end{aligned}
```

In this chain complex framework:

- Standard parity checks are: ``\\partial_3 = H_Z^\\top``, and ``\\partial_2 = H_X``
- Metachecks correspond to: ``\\partial_4 = M_Z^\\top``, and ``\\partial_1 = M_X``

#### Syndrome Validation

The matrices `M_X` and `M_Z` enforce syndrome validity via boundary conditions from the chain complex:

```math
\\begin{aligned}
M_Xs_X = 0 \\quad &\\text{for X-syndromes } (s_X \\in C_1) \\\\
M_Zs_Z = 0 \\quad &\\text{for Z-syndromes } (s_Z \\in C_3)
\\end{aligned}
```

Only CSS codes built using chain complexes and homology have this method.

See also: [`metacheck_matrix_x`](@ref), [`metacheck_matrix_z`](@ref)
"""
function metacheck_matrix end

"""Implemented in a package extension with `Oscar`."""
function bivariate_bicycle_code_k end

"""
    generator_polynomial(c::AbstractCECC)

The generator polynomial g(x) of a [cyclic code](https://en.wikipedia.org/wiki/Cyclic_code)
which generates the ideal corresponding to the code in the quotient ring ``\\mathbb{F}_q[x]/(x^n - 1)``.

The generator polynomial is the unique *monic* polynomial of minimal degree in the
[polynomial code](https://en.wikipedia.org/wiki/Polynomial_code). For a cyclic 
code C of length n over ``\\mathbb{F}_q``, g(x) satisfies:
- g(x) divides ``x^n - 1`` in ``\\mathbb{F}_q[x]``.
- The degree of g(x) is n - k, where k is the code dimension for the non-degenerate case.
- Every codeword polynomial ``c(x) \\in C`` can be expressed as ``c(x) = m(x)g(x) \\mod (x^n - 1)``.

The input is a classical polynomial error-correcting code defined over a finite field.
"""
function generator_polynomial end
