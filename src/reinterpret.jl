import Base: reinterpret

function reinterpret(::Type{U}, p::PauliOperator) where {U<:Unsigned}
    old = eltype(p.xz)
    total_bytes = length(p.xz) * sizeof(old)
    q, r = divrem(total_bytes, sizeof(U))
    if r != 0 || isodd(q)
        throw(ArgumentError(THROW_REINTERPRET_SIZE_MISMATCH))
    end
    new_xz = Base.reinterpret(U, p.xz)
    return PauliOperator(p.phase, p.nqubits, new_xz)
end

function reinterpret(::Type{U}, t::Tableau) where {U<:Unsigned}
    rows = size(t.xzs, 1)
    old = eltype(t.xzs)
    q, r = divrem(rows * sizeof(old), sizeof(U))
    if r != 0 || isodd(q)
        throw(ArgumentError(THROW_REINTERPRET_SIZE_MISMATCH))
    end
    new_xzs = Base.reinterpret(U, t.xzs)
    return Tableau(t.phases, t.nqubits, new_xzs)
end

reinterpret(::Type{U}, s::Stabilizer) where {U<:Unsigned} =
    Stabilizer(reinterpret(U, tab(s)))

reinterpret(::Type{U}, d::Destabilizer) where {U<:Unsigned} =
    Destabilizer(reinterpret(U, tab(d)))

reinterpret(::Type{U}, ms::MixedStabilizer) where {U<:Unsigned} =
    MixedStabilizer(reinterpret(U, tab(ms)), rank(ms))

reinterpret(::Type{U}, md::MixedDestabilizer) where {U<:Unsigned} =
    MixedDestabilizer(reinterpret(U, tab(md)), rank(md))

function reinterpret(::Type{U}, f::PauliFrame) where {U<:Unsigned}
    return PauliFrame(reinterpret(U, f.frame), f.measurements)
end
