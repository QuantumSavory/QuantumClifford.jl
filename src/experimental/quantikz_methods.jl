import Quantikz

function Quantikz.QuantikzOp(op::SparseGate)
    g = op.cliff
    if g==CNOT
        return Quantikz.CNOT(op.indices...)
    elseif g==SWAP*CNOT*SWAP
        return Quantikz.CNOT(op.indices[end:-1:begin]...)
    elseif g==CPHASE
        return Quantikz.CPHASE(op.indices...)
    elseif g==SWAP
        return Quantikz.SWAP(op.indices...)
    else
        return Quantikz.MultiControlU("",[],[],op.indices)
    end
end
Quantikz.QuantikzOp(op::NoisyGate) = Quantikz.QuantikzOp(op.gate)
Quantikz.QuantikzOp(op::BellMeasurement) = Quantikz.ParityMeasurement([string(o) for o in op.pauli], op.indices)
Quantikz.QuantikzOp(op::NoisyBellMeasurement) = Quantikz.QuantikzOp(op.meas)
Quantikz.QuantikzOp(op::VerifyOp) = Quantikz.MultiControlU("\\mathrm{V}",[],[],affectedqubits(op))
Quantikz.QuantikzOp(op::NoiseOp) = Quantikz.Noise(op.indices)
Quantikz.QuantikzOp(op::NoiseOpAll) = Quantikz.NoiseAll()
Quantikz.QuantikzOp(op::Operation) = Quantikz.MultiControlU("",[],[],affectedqubits(op))