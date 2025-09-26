using QuantumClifford:Register, AbstractOperation, applywstatus!, PauliOperator
"""A type of Decision Gate that uses the decoder's lookup table to apply a correction gate based on the syndrome of the input register.
 Does not require an explicit decision function."""
struct DecoderCorrectionGate <: AbstractOperation
    decoderx::AbstractSyndromeDecoder# just a function that maps inputbits to an operation on affectedqubits
    affectedqubits::Vector{Int}
    inputbits::Vector{Int}   
    function DecoderCorrectionGate(dec::TableDecoder, qs, bs)
        return new(dec, collect(Int, qs),collect(Int, bs))
    end
end
function QuantumClifford.apply!(state::Register, op::DecoderCorrectionGate)
  lookup = op.decoder.lookup_table # use the decode function
  key = Vector{bool}(state.bits[op.inputbits])
  if !haskey(lookup, key)
        return state
  end
  n2bit_array = lookup[key]
  pauli_operator = PauliOperator(n2bit_array)
  apply!(state, pauli_operator, op.affectedqubits) 
  state
end
applybranches(s::Register, op::DecoderCorrectionGate; max_order=1) = [(applywstatus!(copy(s),op)...,1,0)]

affectedqubits(op::DecoderCorrectionGate) = op.affected_qubits
affectedbits(op::DecoderCorrectionGate)   = op.input_bits