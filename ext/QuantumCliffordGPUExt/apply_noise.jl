using QuantumClifford: _div, _mod

#according to https://github.com/JuliaGPU/CUDA.jl/blob/ac1bc29a118e7be56d9edb084a4dea4224c1d707/test/core/device/random.jl#L33
#CUDA.jl supports calling rand() inside kernel
function applynoise!(frame::PauliFrameGPU{T},noise::UnbiasedUncorrelatedNoise,i::Int) where {T <: Unsigned}
    p = noise.errprobthird
    lowbit = T(1)
    ibig = _div(T,i-1)+1
    ismall = _mod(T,i-1)
    ismallm = lowbit<<(ismall)

    stab = frame.frame
    xzs = tab(stab).xzs
    rows = size(stab, 1)

    @run_cuda applynoise_kernel(xzs, p, ibig, ismallm, rows) rows
    return frame
end


function applynoise_kernel(xzs::DeviceMatrix{Tme},
    p::Real,
    ibig::Int,
    ismallm::Tme,
    rows::Int) where {Tme <: Unsigned} 

    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end

    r = rand()
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
