import Quantikz

function Quantikz.QuantikzOp(op::SparseGate)
    g = op.cliff
    if g==CNOT
        return Quantikz.CNOT(op.indices...)
    elseif g==SWAP*CNOT*SWAP
        return Quantikz.CNOT(op.indices[end:-1:1]...) # TODO1.3 change to end:-1:begin when you drop julia1.3
    elseif g==CPHASE
        return Quantikz.CPHASE(op.indices...)
    elseif g==SWAP
        return Quantikz.SWAP(op.indices...)
    else
        return Quantikz.MultiControlU("",[],[],op.indices)
    end
end
Quantikz.QuantikzOp(op::AbstractOperation) = Quantikz.MultiControlU(affectedqubits(op))
Quantikz.QuantikzOp(op::BellMeasurement) = Quantikz.ParityMeasurement([string(o) for o in op.pauli], op.indices)
Quantikz.QuantikzOp(op::NoisyBellMeasurement) = Quantikz.QuantikzOp(op.meas)
Quantikz.QuantikzOp(op::ConditionalGate) = Quantikz.ClassicalDecision(affectedqubits(op),op.controlbit)
Quantikz.QuantikzOp(op::DecisionGate) = Quantikz.ClassicalDecision(affectedqubits(op),1) #TODO permit Quantikz.ibegin:Quantikz.iend in ClassicalDecision (and in Noise)
Quantikz.QuantikzOp(op::DenseGate) = Quantikz.MultiControlU(affectedqubits(op))
Quantikz.QuantikzOp(op::DenseMeasurement) = Quantikz.Measurement("\\begin{smallmatrix}$(lstring(op.pauli))\\end{smallmatrix}",affectedqubits(op),op.storagebit)
Quantikz.QuantikzOp(op::NoisyGate) = Quantikz.QuantikzOp(op.gate)
Quantikz.QuantikzOp(op::VerifyOp) = Quantikz.MultiControlU("\\begin{smallmatrix}\\mathrm{Verify:}\\\\$(lstring(op.good_state))\\end{smallmatrix}",affectedqubits(op))
Quantikz.QuantikzOp(op::NoiseOp) = Quantikz.Noise(op.indices)
Quantikz.QuantikzOp(op::NoiseOpAll) = Quantikz.NoiseAll()

function lstring(pauli::PauliOperator)
    v = join(("\\mathtt{$(o)}" for o in replace(string(pauli)[3:end],"_"=>"I")),"\\\\")
end

function lstring(stab::Stabilizer)
    v = join(("\\mathtt{$(replace(string(p),"_"=>"I"))}" for p in stab),"\\\\")
end

Base.show(io::IO, mime::MIME"image/png", circuit::AbstractVector{<:AbstractOperation}; scale=1, kw...) = 
    show(io, mime, [Quantikz.QuantikzOp(c) for c in circuit]; scale=scale, kw...)    
Base.show(io::IO, mime::MIME"image/png", gate::T; scale=1, kw...) where T<:AbstractOperation = 
    show(io, mime, Quantikz.QuantikzOp(gate); scale=scale, kw...)