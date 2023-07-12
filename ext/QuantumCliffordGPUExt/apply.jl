Base.@propagate_inbounds function getxbit(xzs::CuDeviceMatrix{T, 1}, r::Int, c::Int)::T where {T <: Unsigned}
    xzs[QuantumClifford.getbigindex(T, c),r] & QuantumClifford.getmask(T, c)
end
Base.@propagate_inbounds function getzbit(xzs::CuDeviceMatrix{T, 1}, r::Int, c::Int)::T where {T <: Unsigned}
    xzs[end÷2+QuantumClifford.getbigindex(T, c),r]& QuantumClifford.getmask(T, c)
end
Base.@propagate_inbounds function setxbit(xzs::CuDeviceMatrix{T, 1}, r::Int, c::Int, x::T) where {T <: Unsigned}
    cbig = QuantumClifford.getbigindex(T, c)
    xzs[cbig,r] &= ~QuantumClifford.getmask(T, c)
    xzs[cbig,r] |= x
end
Base.@propagate_inbounds function setzbit(xzs::CuDeviceMatrix{T, 1}, r::Int, c::Int, z::T) where {T <: Unsigned}
    cbig = QuantumClifford.getbigindex(T, c)
    xzs[end÷2+cbig,r] &= ~QuantumClifford.getmask(T, c)
    xzs[end÷2+cbig,r] |= z
end
Base.@propagate_inbounds setxbit(xzs::CuDeviceMatrix{T, 1}, r::Int, c::Int, x::T, shift::Int) where {T <: Unsigned} = setxbit(xzs, r, c, x<<shift)
Base.@propagate_inbounds setzbit(xzs::CuDeviceMatrix{T, 1}, r::Int, c::Int, z::T, shift::Int) where {T <: Unsigned} = setzbit(xzs, r, c, z<<shift)

# todo put back the generic types later
# Questions:
# 1- couldn't input tabeulu to gpu kernel
# 2- doesn't support multimodal so I had to write functions one by one
# 3- how to use the getxbit, setxbit functions that are in QuantumClifford? (without having to copy)
# 4- CuArray becomes CuDeviceMatrix in kernel!
function single_qubit_gpu_kernel(xzs::CuDeviceMatrix{Tme, 1},
                                 phases::CuDeviceVector{Tmz, 1},
                                 op::SingleQubitOperator,
                                 rows::Unsigned,
                                 compute_phases::Bool) where {Tme <: Unsigned, Tmz <: Unsigned}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > rows
        return nothing
    end
    c = op.q
    r = idx
    sh = QuantumClifford.getshift(Tme, c)
    xx,zx,xz,zz = Tme.((op.xx,op.zx,op.xz,op.zz)) .<< sh # maybe do this in parent class?
    anticom = ~iszero((~zz & xz & ~xx & zx) | ( zz & ~xz & xx & zx) | (zz &  xz & xx & ~zx))

    # todo. in future each gpu core can be responsible for multiple rows
    x = getxbit(xzs, r, c)
    z = getzbit(xzs, r, c)
    setxbit(xzs, r, c, (x&xx)⊻(z&zx))
    setzbit(xzs, r, c, (x&xz)⊻(z&zz))

    if compute_phases
        if op.px && ~iszero(x)
            phases[r] += 0x2
            phases[r] &= 3
        end
        if op.pz && ~iszero(z)
            phases[r] += 0x2
            phases[r] &= 3
        end
        if ~iszero(x&z) && anticom
            phases[r] += 0x2
            phases[r] &= 3
        end
    end
    return nothing
end

function _apply!(stab::QuantumClifford.Stabilizer{QuantumClifford.Tableau{Tz, Tm}},
    op::QuantumClifford.SingleQubitOperator;
    phases::Val{B}=Val(true)) where {B, Tz<:CuArray{<:Unsigned, 1}, Tm<:CuArray{<:Unsigned, 2}}
    # todo how to use phases similar to before in kernel functions??!
    threads_count = 1024 # Change this later
    rows::Unsigned = size(stab, 1)
    blocks_count = ceil(Int, rows/threads_count)
    tab = QuantumClifford.tab(stab)
    # todo. why can't I pass phases=compute_phases normally without function call?
    CUDA.@sync @cuda threads=threads_count blocks=blocks_count single_qubit_gpu_kernel(tab.xzs, tab.phases, op, rows, B)
end

function two_qubit_gpu_kernel(xzs::CuDeviceMatrix{Tme, 1},
                              phases::CuDeviceVector{Tze, 1},
                              gate::QuantumClifford.AbstractTwoQubitOperator,
                              rows::Unsigned,
                              compute_phases::Bool=true) where {Tme <: Unsigned, Tze <: Unsigned}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > rows
        return nothing
    end

    q1 = gate.q1
    q2 = gate.q2
    shift = QuantumClifford.getshift(Tme, q1) - QuantumClifford.getshift(Tme, q2)
    r = idx;

    _x1::Tme = getxbit(xzs, r, q1)
    _z1::Tme = getzbit(xzs, r, q1)
    _x2::Tme = getxbit(xzs, r, q2)<<shift
    _z2::Tme = getzbit(xzs, r, q2)<<shift
    x1::Tme,z1::Tme,x2::Tme,z2::Tme,phase::Bool = QuantumClifford.qubit_kernel(gate,_x1,_z1,_x2,_z2) # Most `qubit_kernel` functions are defined by a `qubitop2` macro
    setxbit(xzs, r, q1, x1, 0)
    setzbit(xzs, r, q1, z1, 0)
    setxbit(xzs, r, q2, x2, -shift)
    setzbit(xzs, r, q2, z2, -shift)
    if compute_phases && phase
        phases[r] += 0x2
        phases[r] &= 3
    end
    return nothing
end


function _apply!(stab::QuantumClifford.Stabilizer, 
                 gate::G; 
                 phases::Val{B}=Val(true)) where {B, G<:QuantumClifford.AbstractTwoQubitOperator}
    threads_count = 1024 # Change this later
    rows::Unsigned = size(stab, 1)
    blocks_count = ceil(Int, rows/threads_count)
    tab = QuantumClifford.tab(stab)
    # todo. why can't I pass compute_phases=compute_phases normally without function call?
    CUDA.@sync @cuda threads=threads_count blocks=blocks_count two_qubit_gpu_kernel(tab.xzs, tab.phases, gate, rows, B)

    # todo dry this out...!
end
