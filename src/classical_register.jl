"""A register, representing the state of a computer including both a tableaux and an array of classical bits (e.g. for storing measurement results)"""
struct Register{Tzv<:AbstractVector{UInt8},Tm<:AbstractMatrix{<:Unsigned}} <: AbstractQCState
    stab::MixedDestabilizer{Tzv,Tm}
    bits::Vector{Bool}
    Register(s::MixedDestabilizer{A,B}, bits) where {A,B} = new{A,B}(s,bits)
end

Register(s,bits) = Register(MixedDestabilizer(s), bits)

Base.copy(r::Register) = Register(copy(r.stab),copy(r.bits))

stabilizerview(r::Register) = stabilizerview(quantumstate(r))
destabilizerview(r::Register) = destabilizerview(quantumstate(r))
logicalxview(r::Register) = logicalxview(quantumstate(r))
logicalzview(r::Register) = logicalzview(quantumstate(r))

"""A view of the classical bits stored with the state"""
function bitview end
bitview(s::AbstractStabilizer) = ()
bitview(r::Register) = r.bits

"""Only the quantum part of the state (excluding classical bits)"""
function quantumstate end
quantumstate(s::AbstractStabilizer) = s
quantumstate(r::Register) = r.stab

function apply!(r::Register, args...; kwargs...)
    apply!(quantumstate(r), args...; kwargs...)
    r
end
function apply!(r::Register, op::AbstractCliffordOperator, indices; kwargs...)
    apply!(quantumstate(r), op, indices; kwargs...)
    r
end

function apply!(r::Register, m::sMX{T}) where T
    _, res = projectXrand!(quantumstate(r),m.qubit)
    T==Int && (bitview(r)[m.bit] = iszero(res))
    r
end
function apply!(r::Register, m::sMY{T}) where T
    _, res = projectYrand!(quantumstate(r),m.qubit)
    T==Int && (bitview(r)[m.bit] = iszero(res))
    r
end
function apply!(r::Register, m::sMZ{T}) where T
    _, res = projectZrand!(quantumstate(r),m.qubit)
    T==Int && (bitview(r)[m.bit] = iszero(res))
    r
end
function apply!(r::Register, m::PauliMeasurement{A,B,T}) where {A,B,T}
    _, res = projectrand!(quantumstate(r),m.pauli)
    T==Int && (bitview(r)[m.storagebit] = iszero(res))
    r
end
