function apply!(f::QuantumClifford.PauliFrame, op::QuantumClifford.AbstractCliffordOperator)
    _apply!(f.frame, op; phases=Val(false))
    return f
end

function apply_sMZ_kernel!(xzs::CuDeviceMatrix{Tme, 1},
                          measurements::CuDeviceMatrix{Bool, 1},
                          op::sMZ,
                          ibig::Int,
                          ismallm::Tme,
                          rows::Int) where {Tme <: Unsigned} 
    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end    
    should_flip = !iszero(xzs[ibig,f] & ismallm)
    measurements[f,op.bit] = should_flip
    return nothing
end

function apply!(frame::QuantumClifford.PauliFrame, op::QuantumClifford.sMZ) # TODO sMX, sMY
    op.bit == 0 && return frame
    i = op.qubit
    xzs = frame.frame.tab.xzs
    T = eltype(xzs)
    lowbit = T(1)
    ibig = QuantumClifford._div(T,i-1)+1
    ismall = QuantumClifford._mod(T,i-1)
    ismallm = lowbit<<(ismall)


    threads_count = 1024 # todo choose this wiser
    blocks_count = ceil(Int, length(frame)/threads_count)
    # todo make this some kind of macro so that we don't need to write this all the time?!
    CUDA.@sync @cuda threads=threads_count blocks=blocks_count apply_sMZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame))

    return frame
end

function apply_sMRZ_kernel!(xzs::CuDeviceMatrix{Tme, 1},
                          measurements::CuDeviceMatrix{Bool, 1},
                          op::QuantumClifford.sMRZ,
                          ibig::Int, # todo change to Int
                          ismallm::Tme,
                          rows::Int) where {Tme <: Unsigned} 
    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end
    if op.bit != 0 # not good practice. how to replace 0?
        should_flip = !iszero(xzs[ibig,f] & ismallm)
        measurements[f,op.bit] = should_flip
    end
    xzs[ibig,f] &= ~ismallm
    rand(Bool) && (xzs[end÷2+ibig,f] ⊻= ismallm)
    return nothing
end

function apply!(frame::QuantumClifford.PauliFrame, op::QuantumClifford.sMRZ) # TODO sMRX, sMRY
    i = op.qubit
    xzs = frame.frame.tab.xzs
    T = eltype(xzs)
    lowbit = T(1)
    ibig = QuantumClifford._div(T,i-1)+1
    ismall = QuantumClifford._mod(T,i-1)
    ismallm = lowbit<<(ismall)

    threads_count = 1024 # todo choose this wiser
    blocks_count = ceil(Int, length(frame)/threads_count)
    # todo make this some kind of macro so that we don't need to write this all the time?!
    CUDA.@sync @cuda threads=threads_count blocks=blocks_count apply_sMRZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame))
    return frame
end

# todo remove this function later after adding support for CompactifiedCircuit
function pftrajectories(circuit;trajectories=5000)
    ccircuit = if eltype(circuit) <: QuantumClifford.CompactifiedGate
        circuit
    else
        circuit
        ################
        # Apply Compactify after adding apply! support for it...
        # QuantumClifford.compactify_circuit(circuit)
    end

    ccircuit = map(normalize_gate, ccircuit); # todo. remove this. there should be something like this embedded in QuantumClifford

    qmax=maximum((maximum(QuantumClifford.affectedqubits(g)) for g in ccircuit))
    bmax=maximum((maximum(QuantumClifford.affectedbits(g),init=1) for g in ccircuit))
    frames = QuantumClifford.PauliFrame(trajectories, qmax, bmax)
    frames = to_gpu(frames) # todo we can construct this on gpu in the first place. not moving it later...
    pftrajectories(frames, ccircuit)
    return frames
end

function pftrajectories(state::QuantumClifford.PauliFrame, circuit)
    for op in circuit
        apply!(state, op)
    end
    return state
end

# todo this is just a patch. remove it in the future...
normalize_gate(gate::QuantumClifford.AbstractOperation) = gate
normalize_gate(gate::sHadamard) = SingleQubitOperator(gate)
normalize_gate(gate::sPhase) = SingleQubitOperator(gate)
normalize_gate(gate::sInvPhase) = SingleQubitOperator(gate)
normalize_gate(gate::sId1) = SingleQubitOperator(gate)
normalize_gate(gate::sX) = SingleQubitOperator(gate)
normalize_gate(gate::sY) = SingleQubitOperator(gate)
normalize_gate(gate::sZ) = SingleQubitOperator(gate)
