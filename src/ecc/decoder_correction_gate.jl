using QuantumClifford:Register, AbstractOperation, applywstatus!, PauliOperator
import QuantumClifford:affectedqubits,affectedbits,applybranches
"""Apply a Pauli correction determined by a decoder, using syndrome bits from the classical register to correct the specified data qubits."""
struct DecoderCorrectionGate <: AbstractOperation
    decoder::AbstractSyndromeDecoder# just a function that maps inputbits to an operation on affectedqubits
    data_qubits::Vector{Int}
    syndrome_bits::Vector{Int}   
    
    function DecoderCorrectionGate(dec::AbstractSyndromeDecoder, data_qubits, syndrome_bits)
        qs = collect(Int, data_qubits)
        bs = collect(Int, syndrome_bits)
        bits, qubits = size(parity_checks(dec))
        if length(qs)!=qubits
            throw(ArgumentError(lazy"Decoder expects $qubits qubits, got $(length(data_qubits))"))
        end
        if length(bs)!=bits
            throw(ArgumentError(lazy"Decoder expects $bits bits, got $(length(syndrome_bits))"))
        end
        return new(dec, qs,bs)
    end
end

function QuantumClifford.apply!(state::Register, op::DecoderCorrectionGate)
    targets = op.data_qubits
    key = Vector{Bool}(state.bits[op.syndrome_bits])
    all(iszero, key) && return state
    correction = decode(op.decoder, key)  
    correction === nothing && return state
    pauli_operator = PauliOperator(correction)
    apply!(state, pauli_operator, targets) 
    return state
    end

applybranches(s::Register, op::DecoderCorrectionGate; max_order=1) = [(applywstatus!(copy(s),op)...,1,0)]

affectedqubits(op::DecoderCorrectionGate) = op.data_qubits
affectedbits(op::DecoderCorrectionGate)   = op.syndrome_bits

