"""
$(TYPEDEF)

This is a wrapper around a tableau. This "frame" tableau is not to be viewed as a normal stabilizer tableau,
although it does conjugate the same under Clifford operations.
Each row in the tableau refers to a single frame.
The row represents the Pauli operation by which the frame and the reference differ.
"""
struct PauliFrame{T,S} <: AbstractQCState
    frame::T # TODO this should really be a Tableau, but now most of the code seems to be assuming a Stabilizer
    measurements::S # TODO check if when looping over this we are actually looping over the fast axis
end

nqubits(f::PauliFrame) = nqubits(f.frame)
Base.length(f::PauliFrame) = size(f.measurements, 1)
Base.eachindex(f::PauliFrame) = 1:length(f)
Base.copy(f::PauliFrame) = PauliFrame(copy(f.frame), copy(f.measurements))
Base.view(frame::PauliFrame, r) = PauliFrame(view(frame.frame, r), view(frame.measurements, r, :))

fastrow(s::PauliFrame) = PauliFrame(fastrow(s.frame), s.measurements)
fastcolumn(s::PauliFrame) = PauliFrame(fastcolumn(s.frame), s.measurements)

"""
$(TYPEDSIGNATURES)

Prepare an empty set of Pauli frames with the given number of `frames` and `qubits`. Preallocates spaces for `measurement` number of measurements.
"""
function PauliFrame(frames, qubits, measurements)
    stab = fastcolumn(zero(Stabilizer, frames, qubits)) # TODO this should really be a Tableau
    bits = zeros(Bool, frames, measurements)
    frame = PauliFrame(stab, bits)
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
    _apply!(f.frame, op; phases=Val(false))
    return f
end

function apply!(frame::PauliFrame, xor::ClassicalXOR)
    for f in eachindex(frame)
        value = frame.measurements[f,xor.bits[1]]
        for i in xor.bits[2:end]
            value ⊻= frame.measurements[f,i]
        end
        frame.measurements[f, xor.store] = value
    end
end

function apply!(frame::PauliFrame, op::sMX) # TODO implement a faster direct version
    op.bit == 0 && return frame
    apply!(frame, sHadamard(op.qubit))
    apply!(frame, sMZ(op.qubit, op.bit))
end

function apply!(frame::PauliFrame, op::sMRX) # TODO implement a faster direct version
    apply!(frame, sHadamard(op.qubit))
    apply!(frame, sMRZ(op.qubit, op.bit))
end

function apply!(frame::PauliFrame, op::sMZ) # TODO sMY, and faster sMX
    op.bit == 0 && return frame
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

function apply!(frame::PauliFrame, op::sMRZ) # TODO sMRY, and faster sMRX
    i = op.qubit
    xzs = frame.frame.tab.xzs
    T = eltype(xzs)
    lowbit = T(1)
    ibig = _div(T,i-1)+1
    ismall = _mod(T,i-1)
    ismallm = lowbit<<(ismall)

    if op.bit != 0
        @inbounds @simd for f in eachindex(frame)
            should_flip = !iszero(xzs[ibig,f] & ismallm)
            frame.measurements[f,op.bit] = should_flip
        end
    end
    @inbounds @simd for f in eachindex(frame)
        xzs[ibig,f] &= ~ismallm
        rand(Bool) && (xzs[end÷2+ibig,f] ⊻= ismallm)
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

The main method for running Pauli frame simulations of circuits.
See the other methods for lower level access.

Multithreading is enabled by default, but can be disabled by setting `threads=false`.
Do not forget to launch Julia with multiple threads enabled, e.g. `julia -t4`, if you want
to use multithreading.

See also: [`mctrajectories`](@ref), [`petrajectories`](@ref)
"""
function pftrajectories(circuit;trajectories=5000,threads=true)
    _pftrajectories(circuit;trajectories,threads)
end

function _create_pauliframe(ccircuit; trajectories=5000)
    qmax=maximum((maximum(affectedqubits(g)) for g in ccircuit))
    bmax=maximum((maximum(affectedbits(g),init=1) for g in ccircuit))
    return PauliFrame(trajectories, qmax, bmax)
end

function _pftrajectories(circuit;trajectories=5000,threads=true)
    ccircuit = if eltype(circuit) <: CompactifiedGate
        circuit
    else
        compactify_circuit(circuit)
    end
    frames = _create_pauliframe(ccircuit; trajectories)
    nthr = min(Threads.nthreads(),trajectories÷(100))
    if threads && nthr>1
        batchsize = trajectories÷nthr
        Threads.@threads for i in 1:nthr
            b = (i-1)*batchsize+1
            e = i==nthr ? trajectories : i*batchsize
            pftrajectories((@view frames[b:e]), ccircuit)
        end
    else
        pftrajectories(frames, ccircuit)
    end
    return frames
end

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
