using QuantumClifford: get_bitmask_idxs

#according to https://github.com/JuliaGPU/CUDA.jl/blob/ac1bc29a118e7be56d9edb084a4dea4224c1d707/test/core/device/random.jl#L33
#CUDA.jl supports calling rand() inside kernel
function applynoise!(frame::PauliFrameGPU{T},noise::UnbiasedUncorrelatedNoise,i::Int) where {T <: Unsigned}
    xzs = frame.frame.tab.xzs
    p = noise.p
    lowbit, ibig, ismall, ismallm = get_bitmask_idxs(xzs,i)
    rows = size(stab, 1)

    @run_cuda applynoise_kernel(xzs, p, ibig, ismallm, rows) rows
    return frame
end


function applynoise_kernel(xzs::DeviceMatrix{Tₘₑ},
    p::Real,
    ibig::Int,
    ismallm::Tₘₑ,
    rows::Int) where {Tₘₑ <: Unsigned}

    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end

    r = rand()
    p = p/3
    if  r < p # X error
        xzs[ibig,f] ⊻= ismallm
    elseif r < 2p # Z error
        xzs[end÷2+ibig,f] ⊻= ismallm
    elseif r < 3p # Y error
        xzs[ibig,f] ⊻= ismallm
        xzs[end÷2+ibig,f] ⊻= ismallm
    end
    return nothing
end
