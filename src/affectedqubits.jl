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
affectedbits(m::Union{sMRZ,sMZ,sMRX,sMX,sMRY,sMY}) = m.bit==0 ? () : (m.bit,)
affectedbits(c::ClassicalXOR) = (c.bits..., c.store)

function _sentinel_affectedqubits end
QuantumClifford._sentinel_affectedqubits(x::Any) = QuantumClifford.affectedqubits(x)
QuantumClifford._sentinel_affectedqubits(::QuantumClifford.NoiseOpAll) = missing

function _sentinel_maximum end
QuantumClifford._sentinel_maximum(x::Any) = maximum(x)
QuantumClifford._sentinel_maximum(::Missing) = 0