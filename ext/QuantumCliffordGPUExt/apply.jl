using QuantumClifford: getxbit, getzbit, setxbit, setzbit

##############################
# SingleQubitOperator Kernel
##############################

function single_qubit_gpu_kernel(xzs::DeviceMatrix{Tₘₑ},
                                 phases::DeviceVector{Tₚₑ},
                                 op::SingleQubitOperator,
                                 rows::Int,
                                 compute_phases::Bool) where {Tₘₑ <: Unsigned, Tₚₑ <: Unsigned}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > rows
        return nothing
    end
    c = op.q
    r = idx
    sh = QuantumClifford.getshift(Tₘₑ, c)
    xx,zx,xz,zz = Tₘₑ.((op.xx,op.zx,op.xz,op.zz)) .<< sh # maybe putting this in parent function call will improve speed?
    anticom = ~iszero((~zz & xz & ~xx & zx) | ( zz & ~xz & xx & zx) | (zz &  xz & xx & ~zx))

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

##############################
# AbstractSingleQubitOperator Kernel
##############################

function abstract_single_qubit_gpu_kernel(xzs::DeviceMatrix{Tₘₑ},
                                 phases::DeviceVector{Tₚₑ},
                                 gate::QuantumClifford.AbstractSingleQubitOperator,
                                 rows::Int,
                                 compute_phases::Bool) where {Tₘₑ <: Unsigned, Tₚₑ <: Unsigned}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > rows
        return nothing
    end
    c = gate.q
    r = idx

    x::Tₘₑ = getxbit(xzs, r, c)
    z::Tₘₑ = getzbit(xzs, r, c)
    x,z,phase::Bool = QuantumClifford.qubit_kernel(gate,x,z)
    setxbit(xzs, r, c, x)
    setzbit(xzs, r, c, z)
    compute_phases && phase && (phases[r] = (phases[r]+0x2)&3)
    return nothing
end

##############################
# AbstractTwoQubitOperator Kernel
##############################

function two_qubit_gpu_kernel(xzs::DeviceMatrix{Tₘₑ},
                              phases::DeviceVector{Tₚ},
                              gate::QuantumClifford.AbstractTwoQubitOperator,
                              rows::Int,
                              compute_phases::Bool=true) where {Tₘₑ <: Unsigned, Tₚ <: Unsigned}
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > rows
        return nothing
    end

    q1 = gate.q1
    q2 = gate.q2
    shift = QuantumClifford.getshift(Tₘₑ, q1) - QuantumClifford.getshift(Tₘₑ, q2)
    r = idx;

    _x1::Tₘₑ = getxbit(xzs, r, q1)
    _z1::Tₘₑ = getzbit(xzs, r, q1)
    _x2::Tₘₑ = getxbit(xzs, r, q2)<<shift
    _z2::Tₘₑ = getzbit(xzs, r, q2)<<shift
    x1::Tₘₑ,z1::Tₘₑ,x2::Tₘₑ,z2::Tₘₑ,phase::Bool = QuantumClifford.qubit_kernel(gate,_x1,_z1,_x2,_z2) # Most `qubit_kernel` functions are defined by a `qubitop2` macro
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

##############################
# Wrapper _apply! functions
##############################

function _apply!(stab::StabilizerGPU{T},
    op::QuantumClifford.SingleQubitOperator;

    phases::Val{B}=Val(true)) where {B, T <: Unsigned}
    rows = size(stab, 1)
    tab = QuantumClifford.tab(stab)
    (@run_cuda single_qubit_gpu_kernel(tab.xzs, tab.phases, op, rows, B) rows)
    return stab
end

function _apply!(stab::StabilizerGPU{T},
    op::QuantumClifford.AbstractSingleQubitOperator;
    phases::Val{B}=Val(true)) where {B, T <: Unsigned}

    rows = size(stab, 1)
    tab = QuantumClifford.tab(stab)
    (@run_cuda abstract_single_qubit_gpu_kernel(tab.xzs, tab.phases, op, rows, B) rows)
    return stab
end

function _apply!(stab::StabilizerGPU{T},
                 gate::G;
                 phases::Val{B}=Val(true)) where {B, G<:QuantumClifford.AbstractTwoQubitOperator, T <: Unsigned}

    rows = size(stab, 1)
    tab = QuantumClifford.tab(stab)
    (@run_cuda two_qubit_gpu_kernel(tab.xzs, tab.phases, gate, rows, B) rows)
    return stab
end
