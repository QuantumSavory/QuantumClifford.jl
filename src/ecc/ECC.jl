module ECC

using QuantumClifford
using QuantumClifford: AbstractOperation

abstract type AbstractECC end

"""The encoding circuit of a given code."""
function encoding_circuit end

"""Parity check tableau of a code."""
function parity_checks end

"""The number of physical qubits in a code."""
function code_n end

"""The number of stabilizer checks in a code."""
function code_s(c::AbstractECC)
    s = length(parity_checks(c))
    return s
end

"""The number of logical qubits in a code."""
function code_k(c::AbstractECC)
    k = code_n(c) - code_s(c)
    return k
end

"""The rate of a code."""
function rate(c::AbstractECC)
    rate = code_k(c)//code_s(c)
    return rate
end

"""Number of physical qubits for a given parity check tableau"""
function code_n(parity_check_tableau)
    return size(parity_check_tableau)[2]
end

"""Wrapper function for codes of type AbstractECC"""
function naive_syndrome_circuit(code_type::AbstractECC)
    naive_syndrome_circuit(parity_checks(code_type))
end

"""Naive syndrome circuit"""
function naive_syndrome_circuit(parity_check_tableau)
    naive_sc = AbstractOperation[]

    ancilla_bit = 1
    # ancilla_qubit = code_n(c) + ancilla_bit
    for check in parity_check_tableau
        ancilla_qubit = code_n(parity_check_tableau) + ancilla_bit
        for qubit in 1: code_n(parity_check_tableau)
            if check[qubit] == (1,0) # TODO simplify this branch using sXCX and similar gates
                push!(naive_sc, sXCX(qubit, ancilla_qubit))
            elseif check[qubit] == (0,1)
                push!(naive_sc, sCNOT(qubit, ancilla_qubit))
            elseif check[qubit] == (1,1)
                push!(naive_sc, sYCX(qubit, ancilla_qubit))
            end
        end
        mz = sMZ(ancilla_qubit, ancilla_bit)
        push!(naive_sc, mz)
        ancilla_bit +=1
    end

    return naive_sc
end

"""The distance of a code."""
function distance end

"""Parity matrix of a code, given as a stabilizer tableau."""
function parity_matrix(c::AbstractECC)
    paritym = stab_to_gf2(parity_checks(c::AbstractECC))
    return paritym
end

"""Logical X operations of a code."""
function logx_ops(c::AbstractECC)
    MixedDest = MixedDestabilizer(parity_checks(c))
    logicalxview(MixedDest)
end

"""Logical Z operations of a code."""
function logz_ops(c::AbstractECC)
    MixedDest = MixedDestabilizer(parity_checks(c))
    logicalzview(MixedDest)
end

# TODO implement isdegenerate

include("./bitflipcode.jl")
include("./shorcode.jl")
include("./steanecode.jl")

end #module
