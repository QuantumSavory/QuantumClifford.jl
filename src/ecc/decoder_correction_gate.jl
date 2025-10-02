using QuantumClifford:Register, AbstractOperation, applywstatus!, PauliOperator
import QuantumClifford:affectedqubits,affectedbits,applybranches
"""A type of Decision Gate that uses the decoder's lookup table to apply a correction gate based on the syndrome of the input register.
 Does not require an explicit decision function."""
struct DecoderCorrectionGate <: AbstractOperation
    decoder::AbstractSyndromeDecoder# just a function that maps inputbits to an operation on affectedqubits
    affected_qubits::Vector{Int}
    input_bits::Vector{Int}   
    
    function DecoderCorrectionGate(dec::AbstractSyndromeDecoder, affected_qubits, input_bits)
        qs = collect(Int, affected_qubits)
        bs = collect(Int, input_bits)
        if length(qs)!=dec.n
            throw(ArgumentError(lazy"Decoder expects $dec.n qubits, got $(length(affected_qubits))"))
        end
        if length(bs)!=dec.s
            throw(ArgumentError(lazy"Decoder expects $dec.s qubits, got $(length(input_bits))"))
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

