abstract type AbstractECC end
abstract type ClassicalCode end

"""Parity check boolean matrix of a code (only the X entries in the tableau, i.e. the checks for Z errors).

Only CSS codes have this method.

See also: [`parity_checks`](@ref)"""
function parity_checks_x(code::AbstractECC)
    throw(lazy"Codes of type $(typeof(code)) do not have separate X and Z parity checks, either because they are not a CSS code and thus inherently do not have separate checks, or because its separate checks are not yet implemented in this library.")
end

"""Parity check boolean matrix of a code (only the Z entries in the tableau, i.e. the checks for X errors).

Only CSS codes have this method.

See also: [`parity_checks`](@ref)"""
function parity_checks_z(code::AbstractECC)
    throw(lazy"Codes of type $(typeof(code)) do not have separate X and Z parity checks, either because they are not a CSS code and thus inherently do not have separate checks, or because its separate checks are not yet implemented in this library.")
end


"""The number of physical qubits in a code."""
function code_n end

"""The number of stabilizer checks in a code. They might not be all linearly independent, thus `code_s >= code_n-code_k`. For the number of linearly independent checks you can use `LinearAlgebra.rank`."""
function code_s end

"""The distance of a code."""
function distance end

"""The rate of a code."""
function rate(c)
    rate = code_k(c)//code_n(c)
    return rate
end