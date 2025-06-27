# using LoopVectorization
using HostCPUFeatures: pick_vector_width
import SIMD
using GPUArraysCore: AbstractGPUArray

"""Nonvectorized version of `mul_left!` used for unit tests."""
function _mul_left_nonvec!(r::AbstractVector{T}, l::AbstractVector{T}; phases::Bool=true) where T<:Unsigned
    rcnt1, rcnt2 = _mul_ordered_nonvec!(r,l; phases)
    rcnt1 ⊻ (rcnt2 << 1)
end

function _mul_ordered_nonvec!(r::AbstractVector{T}, l::AbstractVector{T}; phases::Bool=true) where T<:Unsigned
    if !phases
        r .⊻= l
        return (0, 0)
    end
    cnt1 = zero(T)
    cnt2 = zero(T)
    len = length(l)÷2
    @inbounds @simd for i in 1:len
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        r[i] = newx1 = x1 ⊻ x2
        r[i+len] = newz1 = z1 ⊻ z2
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        cnt2 ⊻= (cnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        cnt1 ⊻= anti_comm
    end
    rcnt1 = count_ones(cnt1)
    rcnt2 = count_ones(cnt2)
    rcnt1, rcnt2
end

#= # LoopVectorization does not support the necessary bit operations
function mul_ordered_lv!(r::AbstractVector{T}, l::AbstractVector{T}; phases::Val{B}=Val(true)) where {T<:Unsigned, B}
    if !B
        r .⊻= l
        return (0, 0)
    end
    cnt1 = zero(T)
    cnt2 = zero(T)
    len = length(l)÷2
    @turbo for i in 1:len
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        newx1 = x1 ⊻ x2
        r[i] = newx1
        newz1 = z1 ⊻ z2
        r[i+len] = newz1
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        cnt2 ⊻= (newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        cnt1 ⊻= anti_comm
    end
    rcnt1 = count_ones(cnt1)
    rcnt2 = count_ones(cnt2)
    rcnt1, rcnt2
end
=#

function mul_ordered!(r::SubArray{T,1,P,I2,false}, l::AbstractVector{T}; phases::Val{B}=Val(true)) where {T<:Unsigned, B, I2, P}
    # This method exists because SIMD.jl requires fast linear indexing
    # (which is not the case for Adjoint,
    # e.g. when we use `fastcolumn`).
    # The `false` in the SubArray parameters stands for
    # "does not support fast linear indexing".
    # See the other ::SubArray method below as well.
    _mul_ordered_nonvec!(r,l; phases=B)
end

function mul_ordered!(r::SubArray{T,1,P,Tuple{I1, I2},true}, l::AbstractVector{T}; phases::Val{B}=Val(true)) where {T<:Unsigned, B, P, I1<:Any, I2<:AbstractUnitRange}
    # This method exists because SIMD.jl requires fast linear indexing
    # that is NOT strided. See the other ::SubArray method above as well.
    _mul_ordered_nonvec!(r,l; phases=B)
end

function mul_ordered!(r::AbstractVector{T}, l::AbstractVector{T}; phases::Val{B}=Val(true)) where {T<:Unsigned, B}
    if !B
        r .⊻= l
        return (0, 0)
    end
    len = length(l)÷2
    veclen = Int(pick_vector_width(T)) # TODO remove the Int cast
    rcnt1 = 0
    rcnt2 = 0
    packs = len÷veclen
    VT = SIMD.Vec{veclen,T}
    if packs>0
        cnt1 = zero(VT)
        cnt2 = zero(VT)
        lane = SIMD.VecRange{veclen}(0)
        @inbounds for i in 1:veclen:(len-veclen+1)
            # JET-XXX The ::VT should not be necessary, but they help with inference
            x1::VT, x2::VT, z1::VT, z2::VT = l[i+lane], r[i+lane], l[i+len+lane], r[i+len+lane]
            r[i+lane] = newx1 = x1 ⊻ x2
            r[i+len+lane] = newz1 = z1 ⊻ z2
            x1z2 = x1 & z2
            anti_comm = (x2 & z1) ⊻ x1z2
            cnt2 ⊻= (cnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
            cnt1 ⊻= anti_comm
        end
        for i in 1:length(cnt1)
            rcnt1 += count_ones(cnt1[i])
        end
        for i in 1:length(cnt2)
            rcnt2 += count_ones(cnt2[i])
        end
    end
    scnt1 = zero(T)
    scnt2 = zero(T)
    @inbounds for i in (packs*veclen+1):len # same code for the pieces that do not fit in a vector at the end (TODO padding would simplify the code)
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        r[i] = newx1 = x1 ⊻ x2
        r[i+len] = newz1 = z1 ⊻ z2
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        scnt2 ⊻= (scnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        scnt1 ⊻= anti_comm
    end
    rcnt1 += count_ones(scnt1) # anti_comm
    rcnt2 += count_ones(scnt2)
    rcnt1, rcnt2
end

function mul_left!(r::AbstractVector{T}, l::AbstractVector{T}; phases::Val{B}=Val(true))::UInt8 where {T<:Unsigned, B}
    rcnt1, rcnt2 = mul_ordered!(r, l; phases=phases)
    return UInt8((rcnt1 ⊻ (rcnt2<<1))&0x3)
end

function mul_right!(l::AbstractVector{T}, r::AbstractVector{T}; phases::Val{B}=Val(true))::UInt8 where {T<:Unsigned, B}
    rcnt1, rcnt2 = mul_ordered!(l, r; phases=phases)
    return UInt8(((rcnt1 ⊻ (rcnt2<<1)) + rcnt1*2)&0x3) # TODO simplify
end

##############################
# On Pauli Operators
##############################

@inline function mul_left!(r::PauliOperator, l::PauliOperator; phases::Val{B}=Val(true)) where B
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The two Pauli operators should have the same length!")) # TODO skip this when @inbounds is set
    s = mul_left!(r.xz, l.xz, phases=phases)
    B && (r.phase[] = (s+r.phase[]+l.phase[])&0x3)
    r
end

@inline function mul_left!(
    r::PauliOperator{R}, l::PauliOperator{L}; phases::Val{B} = Val(true)
    ) where {R <: AbstractGPUArray, L <: AbstractGPUArray, B}

    nqubits(l) == nqubits(r) || throw(
        DimensionMismatch(
            "The two Pauli operators should have the same length!"
            )
        )
    mul_device!(
        r.xz, r.phase, l.xz, l.phase;
        phases = phases, order_right_left = Val(true)
        )
    return r

end

@inline function mul_right!(l::PauliOperator, r::PauliOperator; phases::Val{B}=Val(true)) where B
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The two Pauli operators should have the same length!")) # TODO skip this when @inbounds is set
    s = mul_right!(l.xz, r.xz, phases=phases)
    B && (l.phase[] = (s+r.phase[]+l.phase[])&0x3)
    l
end

@inline function mul_right!(
    l::PauliOperator{L}, r::PauliOperator{R}; phases::Val{B} = Val(true)
    ) where {L <: AbstractGPUArray, R <: AbstractGPUArray, B}

    nqubits(l) == nqubits(r) || throw(
        DimensionMismatch(
            "The two Pauli operators should have the same length!"
            )
        )
    mul_device!(
        l.xz, l.phase, r.xz, r.phase;
        phases = phases, order_right_left = Val(false)
        )
    return l

end

@inline function mul_left!(r::PauliOperator, l::Tableau, i::Int; phases::Val{B}=Val(true)) where B
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The Pauli operator and the Tableau should have the same length!"))
    s = mul_left!(r.xz, (@view l.xzs[:,i]), phases=phases)
    B && (r.phase[] = (s+r.phase[]+l.phases[i])&0x3)
    r
end

@inline function mul_left!(
    r::PauliOperator{R}, l::Tableau{L}, i::Int; phases::Val{B} = Val(true)
    ) where {R <: AbstractGPUArray, L <: AbstractGPUArray, B}

    nqubits(l) == nqubits(r) || throw(
        DimensionMismatch(
            "The Pauli operator and the Tableau should have the same length!"
            )
        )
    mul_device!(
        r.xz, r.phase, (@view l.xzs[:, i]), (@view l.phases[i]);
        phases = phases, order_right_left = Val(true)
        )
    return r

end

@inline mul_left!(r::PauliOperator, l::Stabilizer, i::Int; phases::Val{B}=Val(true)) where B = mul_left!(r, tab(l), i; phases=phases)

@inline function mul_right!(l::PauliOperator, r::Tableau, i::Int; phases::Val{B}=Val(true)) where B
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The Pauli operator and the Tableau should have the same length!"))
    s = mul_right!(l.xz, (@view r.xzs[:,i]), phases=phases)
    B && (l.phase[] = (s+l.phase[]+r.phases[i])&0x3)
    l
end

@inline function mul_right!(
    l::PauliOperator{L}, r::Tableau{R}, i::Int; phases::Val{B} = Val(true)
    ) where {L <: AbstractGPUArray, R <: AbstractGPUArray, B}

    nqubits(l) == nqubits(r) || throw(
        DimensionMismatch(
            "The Pauli operator and the Tableau should have the same length!"
            )
        )
    mul_device!(
        l.xz, l.phase, (@view r.xzs[:, i]), (@view r.phases[i]);
        phases = phases, order_right_left = Val(false)
        )
    return l

end

@inline mul_right!(l::PauliOperator, r::Stabilizer, i::Int; phases::Val{B}=Val(true)) where B = mul_right!(l, tab(r), i; phases=phases)

##############################
# On Tableaux
##############################

@inline function mul_left!(s::Tableau, m, t::Tableau, i; phases::Val{B}=Val(true)) where B
    nqubits(s) == nqubits(t) || throw(DimensionMismatch("The two Tableaux should have the same length!"))
    extra_phase = mul_left!((@view s.xzs[:,m]), (@view t.xzs[:,i]); phases=phases)
    B && (s.phases[m] = (extra_phase+s.phases[m]+s.phases[i])&0x3)
    s
end

@inline function mul_left!(
    s::Tableau{S}, m, t::Tableau{T}, i; phases::Val{B} = Val(true)
    ) where {S <: AbstractGPUArray, T <: AbstractGPUArray, B}

    nqubits(s) == nqubits(t) || throw(
        DimensionMismatch(
            "The two Tableaux should have the same length!"
            )
        )
    mul_device!(
        (@view s.xzs[:, m]), (@view s.phases[m]),
        (@view t.xzs[:, i]), (@view t.phases[i]);
        phases = phases, order_right_left = Val(true)
        )
    return s

end

@inline function mul_left!(s::Tableau, m, i; phases::Val{B}=Val(true)) where B
    extra_phase = mul_left!((@view s.xzs[:,m]), (@view s.xzs[:,i]); phases=phases)
    B && (s.phases[m] = (extra_phase+s.phases[m]+s.phases[i])&0x3)
    s
end

@inline function mul_left!(
    s::Tableau{S}, m, i; phases::Val{B} = Val(true)
    ) where {S <: AbstractGPUArray, B}

    mul_device!(
        (@view s.xzs[:, m]), (@view s.phases[m]),
        (@view s.xzs[:, i]), (@view s.phases[i]);
        phases = phases, order_right_left = Val(true)
        )
    return s

end

@inline function mul_left!(s::Stabilizer, m, i; phases::Val{B}=Val(true)) where B
    mul_left!(tab(s), m, i; phases)
end

@inline function mul_left!(s::Destabilizer, i, j; phases::Val{B}=Val(true)) where B
    t = tab(s)
    mul_left!(t, j, i; phases=Val(false)) # Indices are flipped to preserve commutation constraints
    n = size(t,1)÷2
    mul_left!(t, i+n, j+n; phases=phases)
end

@inline function mul_left!(s::MixedStabilizer, i, j; phases::Val{B}=Val(true)) where B
    mul_left!(tab(s), i, j; phases=phases)
end

@inline function mul_left!(s::MixedDestabilizer, i, j; phases::Val{B}=Val(true)) where B
    t = tab(s)
    mul_left!(t, j, i; phases=Val(false)) # Indices are flipped to preserve commutation constraints
    n = nqubits(t)
    mul_left!(t, i+n, j+n; phases=phases)
end

@inline function mul_left!(s::Tableau, p::PauliOperator; phases::Val{B}=Val(true)) where B # TODO multithread
    nqubits(s)==nqubits(p) || throw(DimensionMismatch("The Pauli operator and the Tableau should have the same length!"))
    @inbounds @simd for m in eachindex(s)
        extra_phase = mul_left!((@view s.xzs[:,m]), p.xz; phases=phases)
        B && (s.phases[m] = (extra_phase+s.phases[m]+p.phase[])&0x3)
    end
    s
end

@inline function mul_left!(
    s::Tableau{S}, p::PauliOperator{P}; phases::Val{B} = Val(true)
    ) where {S <: AbstractGPUArray, P <: AbstractGPUArray, B}

    nqubits(s) == nqubits(p) || throw(
        DimensionMismatch(
            "The Pauli operator and the Tableau should have the same length!"
            )
        )
    mul_device!(
        s.xzs, s.phases, p.xz, p.phase;
        phases = phases, order_right_left = Val(true)
        )
    return s

end

@inline function mul_left!(s::AbstractStabilizer, p::PauliOperator; phases::Val{B}=Val(true)) where B
    mul_left!(tab(s), p; phases=phases)
    s
end

@inline function mul_right!(s::Tableau, p::PauliOperator; phases::Val{B}=Val(true)) where B # TODO multithread
    nqubits(s)==nqubits(p) || throw(DimensionMismatch("The Pauli operator and the Tableau should have the same length!"))
    @inbounds @simd for m in eachindex(s)
        extra_phase = mul_right!((@view s.xzs[:,m]), p.xz; phases=phases)
        B && (s.phases[m] = (extra_phase+s.phases[m]+p.phase[])&0x3)
    end
    s
end

@inline function mul_right!(
    s::Tableau{S}, p::PauliOperator{P}; phases::Val{B} = Val(true)
    ) where {S <: AbstractGPUArray, P <: AbstractGPUArray, B}

    nqubits(s) == nqubits(p) || throw(
        DimensionMismatch(
            "The Pauli operator and the Tableau should have the same length!"
            )
        )
    mul_device!(
        s.xzs, s.phases, p.xz, p.phase;
        phases = phases, order_right_left = Val(false)
        )
    return s

end

@inline function mul_right!(s::AbstractStabilizer, p::PauliOperator; phases::Val{B}=Val(true)) where B
    mul_right!(tab(s), p; phases=phases)
    s
end
