function apply!(f::PauliFrameGPU{T}, op::QuantumClifford.AbstractCliffordOperator) where {T <: Unsigned}
    _apply!(f.frame, op; phases=Val(false))
    return f
end

function apply_sMZ_kernel!(xzs::DeviceMatrix{Tme},
                          measurements::DeviceMatrix{Bool},
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

function apply!(frame::PauliFrameGPU{T}, op::QuantumClifford.sMZ) where {T <: Unsigned} # TODO sMX, sMY
    op.bit == 0 && return frame
    i = op.qubit
    xzs = frame.frame.tab.xzs
    lowbit = T(1)
    ibig = QuantumClifford._div(T,i-1)+1
    ismall = QuantumClifford._mod(T,i-1)
    ismallm = lowbit<<(ismall)
    CUDA.@sync @run_cuda apply_sMZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame)) length(frame)
    return frame
end

function apply_sMRZ_kernel!(xzs::DeviceMatrix{Tme},
                          measurements::DeviceMatrix{Bool},
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

function apply!(frame::PauliFrameGPU{T}, op::QuantumClifford.sMRZ) where {T <: Unsigned} # TODO sMRX, sMRY
    i = op.qubit
    xzs = frame.frame.tab.xzs
    lowbit = T(1)
    ibig = QuantumClifford._div(T,i-1)+1
    ismall = QuantumClifford._mod(T,i-1)
    ismallm = lowbit<<(ismall)
    CUDA.@sync @run_cuda apply_sMRZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame)) length(frame)
    return frame
end
