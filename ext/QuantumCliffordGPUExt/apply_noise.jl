using QuantumClifford: _div, _mod

function applynoise!(frame::PauliFrameGPU{T},noise::UnbiasedUncorrelatedNoise,i::Int) where {T <: Unsigned}
    p = noise.errprobthird
    lowbit = T(1)
    ibig = _div(T,i-1)+1
    ismall = _mod(T,i-1)
    ismallm = lowbit<<(ismall)

    stab = frame.frame
    xzs = tab(stab).xzs
    rows::Unsigned = size(stab, 1)

    @run_cuda applynoise_kernel(xzs, p, ibig, ismallm, rows) rows
    return frame
end


function applynoise_kernel(xzs::DeviceMatrix{Tme},
    p::Real,
    ibig::Int, # todo why is this Int and others Tme?
    ismallm::Tme, # why is rows Unsigned sometimes and sometimes Int?
    rows::Unsigned) where {Tme <: Unsigned} 

    f = (blockIdx().x - 1) * blockDim().x + threadIdx().x;
    if f > rows
        return nothing
    end

    r = rand() # todo add type to it? is it uniform?
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
