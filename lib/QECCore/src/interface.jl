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

The parity check matrix of a error correction code. The matrix is shaped as (X|Z).

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

See also: [`parity_matrix`](@ref)
"""
function parity_matrix_x end

"""
    parity_matrix_z(c::AbstractCSSCode)

Parity check boolean matrix of a code (only the Z entries in the tableau, i.e. the checks for X errors).
Only CSS codes have this method.

See also: [`parity_matrix`](@ref)
"""
function parity_matrix_z end


"""
    code_n(c::AbstractECC)

The number of physical qubits in a error correction code.

See also: [`code_k`](@ref) and [`code_s`](@ref)
"""
function code_n end

"""
    code_s(c::AbstractECC)

The number of stabilizers in a error correction code. They might not be all linearly independent, thus `code_s >= code_n-code_k`. For the number of linearly independent checks you can use `LinearAlgebra.rank`.

See also: [`code_n`](@ref) and [`code_k`](@ref)
"""
function code_s end

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