using QuantumClifford: get_bitmask_idxs

##############################
# sMZ
##############################

function apply_sMZ_kernel!(xzs::DeviceMatrix{Tₘₑ},
                          measurements::DeviceMatrix{Bool},
                          op::sMZ,
                          ibig::Int,
                          ismallm::Tₘₑ,
                          rows::Int) where {Tₘₑ <: Unsigned}
    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end
    should_flip = !iszero(xzs[ibig,f] & ismallm)
    measurements[f,op.bit] = should_flip
    return nothing
end

function apply!(frame::PauliFrameGPU{T}, op::QuantumClifford.sMZ) where {T <: Unsigned}
    op.bit == 0 && return frame
    i = op.qubit
    xzs = frame.frame.tab.xzs
    lowbit, ibig, ismall, ismallm = get_bitmask_idxs(xzs,i)
    (@run_cuda apply_sMZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame)) length(frame))
    return frame
end

##############################
# sMRZ
##############################

function apply_sMRZ_kernel!(xzs::DeviceMatrix{Tₘₑ},
                          measurements::DeviceMatrix{Bool},
                          op::QuantumClifford.sMRZ,
                          ibig::Int, # todo change to Int
                          ismallm::Tₘₑ,
                          rows::Int) where {Tₘₑ <: Unsigned}
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
    lowbit, ibig, ismall, ismallm = get_bitmask_idxs(xzs,i)
    (@run_cuda apply_sMRZ_kernel!(xzs, frame.measurements, op, ibig, ismallm, length(frame)) length(frame))
    return frame
end
