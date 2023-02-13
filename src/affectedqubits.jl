"""A method giving the qubits acted upon by a given operation. Part of the Noise interface."""
function affectedqubits end
affectedqubits(g::AbstractSingleQubitOperator) = [g.q,]
affectedqubits(g::AbstractTwoQubitOperator) = [g.q1, g.q2]
affectedqubits(g::NoisyGate) = affectedqubits(g.gate)
affectedqubits(g::SparseGate) = g.indices
affectedqubits(b::BellMeasurement) = [m.qubit for m in b.measurements]
affectedqubits(r::Reset) = r.indices
affectedqubits(n::NoiseOp) = n.indices
affectedqubits(g::PauliMeasurement) = 1:length(g.pauli)
affectedqubits(m::AbstractMeasurement) = [m.qubit]
