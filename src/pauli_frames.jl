"""
$(TYPEDEF)

The frame field holds a tableau. This tableau is not to be viewed as a normal stabilizer tableau, although
it does conjugate the same under Clifford operations. Each row in the tableau refers to a single frame. Each column a qubit.
The values stored are X,Y,Z,_ and indicate whether the corresponding pauli operation is waiting to be applied to that qubit.

The provided reference measurements will be used to create the measurements field of this struct. They will also be
used when performing measurement. The indices of the provided measurements must match the secondary value of future measurements
when doing sMZ(measure_this_qubit, index_of_reference_provided). See also [`apply!(frame::PauliFrame, op:sMZ)`](@ref)
"""
struct PauliFrame{T} <: AbstractQCState
    frame::T # TODO this should really be a Tableau
    measurements::Matrix{Bool}
end

nqubits(f::PauliFrame) = nqubits(f.frame)
Base.length(f::PauliFrame) = size(f.measurements, 1)
Base.eachindex(f::PauliFrame) = 1:length(f)

"""
$(TYPEDSIGNATURES)

Prepare an empty set of Pauli frames with the given number of `frames` and `qubits`. Preallocates spaces for `measurement` number of measurements.
"""
function PauliFrame(frames, qubits, measurements)
    stab = zero(Stabilizer, frames, qubits) # TODO this should really be a Tableau
    frame = PauliFrame(stab, zeros(Bool, frames, measurements))
    initZ!(frame)
    return frame
end

"""
$(TYPEDSIGNATURES)

Inject random Z errors over all frames and qubits for the supplied PauliFrame with probability 0.5.

Calling this after initialization is essential for simulating any non-deterministic circuit.
It is done automatically by most [`PauliFrame`](@ref) constructors.
"""
function initZ!(frame::PauliFrame)
    T = eltype(frame.frame.tab.xzs)

    @inbounds @simd for f in eachindex(frame)
        @simd for row in 1:size(frame.frame.tab.xzs,1)÷2
            frame.frame.tab.xzs[end÷2+row,f] = rand(T)
        end
    end
    return frame
end

function apply!(f::PauliFrame, op::AbstractCliffordOperator)
    QuantumClifford._apply!(f.frame, op; phases=Val(false))
    return f
end

function apply!(frame::PauliFrame, op::sMZ)
    i = op.qubit
    xzs = frame.frame.tab.xzs
    T = eltype(xzs)
    lowbit = T(1)
    ibig = _div(T,i-1)+1
    ismall = _mod(T,i-1)
    ismallm = lowbit<<(ismall)

    @inbounds @simd for f in eachindex(frame)
        should_flip = !iszero(xzs[ibig,f] & ismallm)
        frame.measurements[f,op.bit] = should_flip
    end

    return frame
end

function applynoise!(frame::PauliFrame,noise::UnbiasedUncorrelatedNoise,i::Int)
    p = noise.errprobthird
    T = eltype(frame.frame.tab.xzs)

    lowbit = T(1)
    ibig = _div(T,i-1)+1
    ismall = _mod(T,i-1)
    ismallm = lowbit<<(ismall)

    @inbounds @simd for f in eachindex(frame)
        r = rand()
        if  r < p # X error
            frame.frame.tab.xzs[ibig,f] ⊻= ismallm
        elseif r < 2p # Z error
            frame.frame.tab.xzs[end÷2+ibig,f] ⊻= ismallm
        elseif r < 3p # Y error
            frame.frame.tab.xzs[ibig,f] ⊻= ismallm
            frame.frame.tab.xzs[end÷2+ibig,f] ⊻= ismallm
        end
    end
    return frame
end

"""
Perform a "Pauli frame" style simulation of a quantum circuit.
"""
function pftrajectories end

"""
$(TYPEDSIGNATURES)

Evolve each frame stored in [`PauliFrame`](@ref) by the given circuit.
"""
function pftrajectories(state::PauliFrame, circuit)
    for op in circuit
        apply!(state, op)
    end
    return state
end

"""
$(TYPEDSIGNATURES)

For a given [`Register`](@ref) and circuit, simulates the reference circuit acting on the
register and then also simulate numerous [`PauliFrame`](@ref) trajectories.
Returns the register and the [`PauliFrame`](@ref) instance.

Use [`pfmeasurements`](@ref) to get the measurement results.
"""
function pftrajectories(register::Register, circuit; trajectories=500)
    for op in circuit
        apply!(register, op)
    end
    frame = PauliFrame(trajectories, nqubits(register), length(bitview(register)))
    pftrajectories(frame, circuit)
    register, frame
end

"""
For a given simulated state, e.g. a [`PauliFrame`](@ref) instance, returns the measurement results.
"""
function pfmeasurement end

"""
$(TYPEDSIGNATURES)

Returns the measurements stored in the bits of the given [`Register`](@ref).
"""
pfmeasurements(register::Register) = bitview(register)

"""
$(TYPEDSIGNATURES)

Returns the measurement results for each frame in the [`PauliFrame`](@ref) instance.

!!! warning "Relative mesurements"
    The return measurements are relative to the reference measurements, i.e. they only say
    whether the reference measurements have been flipped in the given frame.
"""
pfmeasurements(frame::PauliFrame) = frame.measurements

"""
$(TYPEDSIGNATURES)

Takes the references measurements from the given [`Register`](@ref) and applies the flips
as prescribed by the [`PauliFrame`](@ref) relative measurements. The result is the actual
(non-relative) measurement results for each frame.
"""
pfmeasurements(register::Register, frame::PauliFrame) = pfmeasurements(register) .⊻ pfmeasurements(frame)
