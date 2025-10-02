using QuantumClifford:Register, AbstractOperation, applywstatus!, PauliOperator
import QuantumClifford:affectedqubits,affectedbits,applybranches
"""A gate that uses a lookup table to apply a correction gate based on the syndrome of the input register."""
struct DecoderCorrectionGate <: AbstractOperation
    decoder::AbstractSyndromeDecoder# just a function that maps inputbits to an operation on affectedqubits
    affected_qubits::Vector{Int}
    input_bits::Vector{Int}   
    
    function DecoderCorrectionGate(dec::AbstractSyndromeDecoder, affected_qubits, input_bits)
        qs = collect(Int, affected_qubits)
        bs = collect(Int, input_bits)
        bits, qubits = size(parity_checks(decoder))
        if length(qs)!=qubits
            throw(ArgumentError(lazy"Decoder expects $qubits qubits, got $(length(affected_qubits))"))
        end
        if length(bs)!=bits
            throw(ArgumentError(lazy"Decoder expects $bits bits, got $(length(input_bits))"))
        end
        return new(dec, qs,bs)
    end
end
function QuantumClifford.apply!(state::Register, op::DecoderCorrectionGate)
  targets = op.affected_qubits
  key = Vector{Bool}(state.bits[op.input_bits])
  correction = decode(op.decoder, key)  
  correction === nothing && return state
  pauli_operator = PauliOperator(correction)
  apply!(state, pauli_operator, targets) 
  return state
end
applybranches(s::Register, op::DecoderCorrectionGate; max_order=1) = [(applywstatus!(copy(s),op)...,1,0)]

affectedqubits(op::DecoderCorrectionGate) = op.affected_qubits
affectedbits(op::DecoderCorrectionGate)   = op.input_bits

