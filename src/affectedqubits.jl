"""A method giving the qubits acted upon by a given operation. Part of the Noise interface."""
function affectedqubits end
affectedqubits(g::AbstractSingleQubitOperator) = (g.q,)
affectedqubits(g::AbstractTwoQubitOperator) = (g.q1, g.q2)
affectedqubits(g::NoisyGate) = affectedqubits(g.gate,)
affectedqubits(g::SparseGate) = g.indices
affectedqubits(b::BellMeasurement) = map(m->m.qubit, b.measurements)
affectedqubits(r::Reset) = r.indices
affectedqubits(n::NoiseOp) = n.indices
affectedqubits(g::PauliMeasurement) = 1:length(g.pauli)
affectedqubits(p::PauliOperator) = 1:length(p)
affectedqubits(m::Union{AbstractMeasurement,sMRX,sMRY,sMRZ}) = (m.qubit,)
affectedqubits(v::VerifyOp) = v.indices
affectedqubits(c::CliffordOperator) = 1:nqubits(c)
affectedqubits(c::ClassicalXOR) = ()

affectedbits(o) = ()
affectedbits(m::sMRZ) = (m.bit,)
affectedbits(m::sMZ) = (m.bit,)
affectedbits(c::ClassicalXOR) = (c.bits..., c.store)
