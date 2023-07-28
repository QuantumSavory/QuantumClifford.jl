using QuantumClifford: getxbit, getzbit, setxbit, setzbit

function single_qubit_gpu_kernel(xzs::DeviceMatrix{Tme},
                                 phases::DeviceVector{Tmz},
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

function abstract_single_qubit_gpu_kernel(xzs::DeviceMatrix{Tme},
                                 phases::DeviceVector{Tmz},
                                 gate::QuantumClifford.AbstractSingleQubitOperator,
                                 rows::Unsigned,
                                 compute_phases::Bool) where {Tme <: Unsigned, Tmz <: Unsigned}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > rows
        return nothing
    end
    c = gate.q
    r = idx

    x::Tme = getxbit(xzs, r, c)
    z::Tme = getzbit(xzs, r, c)
    x,z,phase::Bool = QuantumClifford.qubit_kernel(gate,x,z)
    setxbit(xzs, r, c, x)
    setzbit(xzs, r, c, z)
    compute_phases && phase && (phases[r] = (phases[r]+0x2)&3)
    return nothing
end

function _apply!(stab::StabilizerGPU{T},
    op::QuantumClifford.SingleQubitOperator;
    phases::Val{B}=Val(true)) where {B, T <: Unsigned}
    # todo how to use phases similar to before in kernel functions??!
    rows::Unsigned = size(stab, 1)
    tab = QuantumClifford.tab(stab)
    # todo. why can't I pass phases=compute_phases normally without function call?
    CUDA.@sync (@run_cuda single_qubit_gpu_kernel(tab.xzs, tab.phases, op, rows, B) rows)
    return stab
end

function _apply!(stab::StabilizerGPU{T},
    op::QuantumClifford.AbstractSingleQubitOperator;
    phases::Val{B}=Val(true)) where {B, T <: Unsigned}

    rows::Unsigned = size(stab, 1)
    tab = QuantumClifford.tab(stab)
    # todo. why can't I pass phases=compute_phases normally without function call?
    CUDA.@sync (@run_cuda abstract_single_qubit_gpu_kernel(tab.xzs, tab.phases, op, rows, B) rows)
    return stab
end

function two_qubit_gpu_kernel(xzs::DeviceMatrix{Tme},
                              phases::DeviceVector{Tze},
                              gate::QuantumClifford.AbstractTwoQubitOperator, # todo. change to two qubit operator instead os abstract version
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


function _apply!(stab::StabilizerGPU{T}, 
                 gate::G; 
                 phases::Val{B}=Val(true)) where {B, G<:QuantumClifford.AbstractTwoQubitOperator, T <: Unsigned} # todo. change to two qubit operator instead os abstract version
    rows::Unsigned = size(stab, 1)
    tab = QuantumClifford.tab(stab)
    # todo. why can't I pass compute_phases=compute_phases normally without function call?
    CUDA.@sync (@run_cuda two_qubit_gpu_kernel(tab.xzs, tab.phases, gate, rows, B) rows)
    return stab
end
