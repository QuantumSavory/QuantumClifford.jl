module QuantumCliffordQuantikzExt

import Quantikz
using QuantumClifford
using QuantumClifford.Experimental.NoisyCircuits
using QuantumClifford: AbstractOperation

function Quantikz.QuantikzOp(op::SparseGate)
    g = op.cliff
    if g==tCNOT
        return Quantikz.CNOT(op.indices...)
    elseif g==tSWAP*tCNOT*tSWAP
        return Quantikz.CNOT(op.indices[end:-1:begin]...)
    elseif g==tCPHASE
        return Quantikz.CPHASE(op.indices...)
    elseif g==tSWAP
        return Quantikz.SWAP(op.indices...)
    else
        return Quantikz.MultiControlU([],[],op.indices) # TODO Permit skipping the string
    end
end
Quantikz.QuantikzOp(op::AbstractOperation) = Quantikz.MultiControlU(affectedqubits(op))
Quantikz.QuantikzOp(op::PauliOperator) = Quantikz.MultiControlU("\\begin{array}{c}$(lstring(op))\\end{array}", affectedqubits(op))
Quantikz.QuantikzOp(op::sId1) = Quantikz.Id(affectedqubits(op)...)
function Quantikz.QuantikzOp(op::AbstractSingleQubitOperator)
    T = typeof(op)
    str = if T in (sHadamard, sPhase, sX, sY, sZ)
        string(T)[2:2]
    elseif T == sInvPhase
        "P^{-1}"
    else
        "U"
    end
    str = "{"*str*"}"
    return Quantikz.U(str,affectedqubits(op)...)
end
Quantikz.QuantikzOp(op::sCNOT) = Quantikz.CNOT(affectedqubits(op)...)
Quantikz.QuantikzOp(op::sZCX) = Quantikz.CNOT(affectedqubits(op)...)
Quantikz.QuantikzOp(op::sXCZ) = Quantikz.CNOT(reverse(affectedqubits(op))...)
Quantikz.QuantikzOp(op::sCPHASE) = Quantikz.CPHASE(affectedqubits(op)...)
Quantikz.QuantikzOp(op::sSWAP) = Quantikz.SWAP(affectedqubits(op)...)
Quantikz.QuantikzOp(op::sZCZ) = Quantikz.CPHASE(affectedqubits(op)...)
Quantikz.QuantikzOp(op::sXCX) = Quantikz.MultiControl([],[],affectedqubits(op),[])
Quantikz.QuantikzOp(op::sZCY) = Quantikz.MultiControlU("Y",affectedqubits(op)...)
Quantikz.QuantikzOp(op::sYCZ) = Quantikz.MultiControlU("Y",reverse(affectedqubits(op))...)
function Quantikz.QuantikzOp(op::AbstractTwoQubitOperator)
    T = typeof(op)
    if T in (sXCY,sYCX,sYCY,sYCZ,sZCY)
        return Quantikz.U(string(T),affectedqubits(op))
    else
        return Quantikz.U("{U}",affectedqubits(op)...)
    end
end
Quantikz.QuantikzOp(op::BellMeasurement) = Quantikz.ParityMeasurement(["\\mathtt{$(string(typeof(o))[3])}" for o in op.measurements], affectedqubits(op))
Quantikz.QuantikzOp(op::NoisyBellMeasurement) = Quantikz.QuantikzOp(op.meas)
Quantikz.QuantikzOp(op::ConditionalGate) = Quantikz.ClassicalDecision(affectedqubits(op),op.controlbit)
Quantikz.QuantikzOp(op::DecisionGate) = Quantikz.ClassicalDecision(affectedqubits(op),Quantikz.ibegin:Quantikz.iend)
#Quantikz.QuantikzOp(op::DenseGate) = Quantikz.MultiControlU(affectedqubits(op))
Quantikz.QuantikzOp(op::PauliMeasurement) = Quantikz.Measurement("\\begin{array}{c}$(lstring(op.pauli))\\end{array}",affectedqubits(op),op.bit)
Quantikz.QuantikzOp(op::sMX) = Quantikz.Measurement("\\begin{array}{c}\\mathtt{X}\\end{array}",affectedqubits(op), iszero(op.bit) ? nothing : op.bit)
Quantikz.QuantikzOp(op::sMY) = Quantikz.Measurement("\\begin{array}{c}\\mathtt{Y}\\end{array}",affectedqubits(op), iszero(op.bit) ? nothing : op.bit)
Quantikz.QuantikzOp(op::sMZ) = Quantikz.Measurement("\\begin{array}{c}\\mathtt{Z}\\end{array}",affectedqubits(op), iszero(op.bit) ? nothing : op.bit)
#Quantikz.QuantikzOp(op::SparseMeasurement) = Quantikz.Measurement("\\begin{array}{c}$(lstring(op.pauli))\\end{array}",affectedqubits(op),op.bit)
Quantikz.QuantikzOp(op::NoisyGate) = Quantikz.QuantikzOp(op.gate)
Quantikz.QuantikzOp(op::VerifyOp) = Quantikz.MultiControlU("\\begin{array}{c}\\mathrm{Verify:}\\\\$(lstring(op.good_state))\\end{array}",affectedqubits(op))
function Quantikz.QuantikzOp(op::Reset) # TODO This is complicated because quantikz requires $$ for some operators but not all of them... Fix in Quantikz.jl
    m,M = extrema(op.indices)
    indices = sort(op.indices)
    str = "\\begin{array}{c}\\\\$(lstring(op.resetto))\\end{array}"
    if collect(m:M)==indices
        Quantikz.Initialize("\$$str\$",affectedqubits(op))
    else
        Quantikz.Initialize("$str",affectedqubits(op))
    end
end
Quantikz.QuantikzOp(op::NoiseOp) = Quantikz.Noise(op.indices)
Quantikz.QuantikzOp(op::NoiseOpAll) = Quantikz.NoiseAll()

function lstring(pauli::PauliOperator)
    v = join(("\\mathtt{$(o)}" for o in replace(string(pauli)[3:end],"_"=>"I")),"\\\\")
end

function lstring(stab::Stabilizer)
    v = join(("\\mathtt{$(replace(string(p),"_"=>"I"))}" for p in stab),"\\\\")
end

function Base.show(io::IO, mime::MIME"image/png", circuit::AbstractVector{<:AbstractOperation}; scale=1, kw...)
    if length(circuit)==0
        show(io, mime, [Quantikz.Id(1)]; scale=scale, kw...)
    else
        show(io, mime, [Quantikz.QuantikzOp(c) for c in circuit]; scale=scale, kw...)
    end
end
Base.show(io::IO, mime::MIME"image/png", gate::T; scale=1, kw...) where T<:AbstractOperation =
    show(io, mime, Quantikz.QuantikzOp(gate); scale=scale, kw...)

end
