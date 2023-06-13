"""A register, representing the state of a computer including both a tableaux and an array of classical bits (e.g. for storing measurement results)"""
struct Register{T<:Tableau} <: AbstractQCState # TODO simplify type parameters (remove nesting)
    stab::MixedDestabilizer{T}
    bits::Vector{Bool}
    Register(s::MixedDestabilizer{T}, bits) where {T} = new{T}(s,bits)
end

Register(s,bits) = Register(MixedDestabilizer(s), bits)
Register(s) = Register(s, Bool[])
Register(s::MixedDestabilizer,nbits::Int) = Register(s, falses(nbits))

Base.copy(r::Register) = Register(copy(r.stab),copy(r.bits))
Base.:(==)(l::Register,r::Register) = l.stab==r.stab && l.bits==r.bits

stabilizerview(r::Register) = stabilizerview(quantumstate(r))
destabilizerview(r::Register) = destabilizerview(quantumstate(r))
logicalxview(r::Register) = logicalxview(quantumstate(r))
logicalzview(r::Register) = logicalzview(quantumstate(r))

nqubits(r::Register) = nqubits(r.stab)

"""A view of the classical bits stored with the state"""
function bitview end
bitview(s::AbstractStabilizer) = ()
bitview(r::Register) = r.bits

"""Only the quantum part of the state (excluding classical bits)"""
function quantumstate end
quantumstate(s::AbstractStabilizer) = s
quantumstate(r::Register) = r.stab

tab(r::Register) = tab(quantumstate(r))

tensor(regs::Register...) = Register(tensor((quantumstate(r) for r in regs)...), [bit for r in regs for bit in r.bits])

function apply!(r::Register, op, args...; kwargs...)
    apply!(quantumstate(r), op, args...; kwargs...)
    r
end

function apply!(r::Register, m::sMX)
    _, res = projectXrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end
function apply!(r::Register, m::sMY)
    _, res = projectYrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end
function apply!(r::Register, m::sMZ)
    _, res = projectZrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end
function apply!(r::Register, m::sMRX) # TODO sMRY
    _, anticom, res = projectX!(quantumstate(r),m.qubit)
    mres = if isnothing(res)
        mres = rand((0x0, 0x2))
        phases(stabilizerview(r))[anticom] = mres
        mres
    else
        res
    end
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(mres))
    if mres==0x2
        apply!(r, sZ(m.qubit))
    end
    r
end
function apply!(r::Register, m::sMRZ) # TODO sMRY
    _, anticom, res = projectZ!(quantumstate(r),m.qubit)
    mres = if isnothing(res)
        mres = rand(Bool)
        phases(stabilizerview(r))[anticom] = 0x2*mres
        mres
    else
        !iszero(res)
    end
    m.bit!=0 && (bitview(r)[m.bit] = mres)
    if mres
        apply!(r, sX(m.qubit))
    end
    r
end
function apply!(r::Register, m::PauliMeasurement{A,B}) where {A,B}
    _, res = projectrand!(r,m.pauli)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end
function projectXrand!(r::Register, m)
    q = quantumstate(r)
    _, res = projectXrand!(q,m)
    r, res
end
function projectYrand!(r::Register, m)
    q = quantumstate(r)
    _, res = projectYrand!(q,m)
    r, res
end
function projectZrand!(r::Register, m)
    q = quantumstate(r)
    _, res = projectZrand!(q,m)
    r, res
end
function projectrand!(r::Register, m)
    q = quantumstate(r)
    _, res = projectrand!(q,m)
    r, res
end

function traceout!(r::Register, arg)
    q = quantumstate(r)
    traceout!(q,arg)
    q
end
