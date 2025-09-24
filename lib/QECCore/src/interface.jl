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

Returns the `X`-metacheck matrix (``\\partial_{Mx}`` boundary map in
[chain complex](https://en.wikipedia.org/wiki/Chain_complex) notation) for a CSS code.

This matrix verifies validity of `X`-syndromes (`Z`-error measurements).

Only CSS codes built using chain complexes and homology have this method.

See also: [`metacheck_matrix_z`](@ref), [`metacheck_matrix`](@ref), [`parity_matrix_x`](@ref)
"""
function metacheck_matrix_x end

"""
    metacheck_matrix_z(c::AbstractCSSCode)

Returns the `Z`-metacheck matrix (``\\partial_{Mz}`` boundary map in
[chain complex](https://en.wikipedia.org/wiki/Chain_complex) notation) for a CSS code.

This matrix verifies validity of `Z`-syndromes (`X`-error measurements).

Only CSS codes built using chain complexes and homology have this method.

See also: [`metacheck_matrix_x`](@ref), [`metacheck_matrix`](@ref), [`parity_matrix_z`](@ref)
"""
function metacheck_matrix_z end

"""
    metacheck_matrix(c::AbstractCSSCode)

Returns both `X` and `Z` metacheck matrices.

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
