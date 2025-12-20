import Base: reinterpret

function reinterpret(::Type{U}, p::PauliOperator) where {U<:Unsigned}
    old = eltype(p.xz)
    total_bytes = length(p.xz) * sizeof(old)
    if total_bytes % sizeof(U) != 0
        throw(
            ArgumentError(
                "cannot reinterpret PauliOperator: backing bytes not divisible by sizeof($(U))",
            ),
        )
    end
    new_xz = Base.reinterpret(U, p.xz)
    nl = length(new_xz)
    if nl % 2 != 0
        throw(
            ArgumentError(
                "cannot reinterpret PauliOperator: resulting backing array length $(nl) is not even (cannot split into X and Z halves).",
            ),
        )
    end
    return PauliOperator(p.phase, p.nqubits, new_xz)
end

function reinterpret(::Type{U}, t::Tableau) where {U<:Unsigned}
    rows = size(t.xzs, 1)
    cols = size(t.xzs, 2)
    old = eltype(t.xzs)
    total_bytes = rows * sizeof(old)
    if total_bytes % sizeof(U) != 0
        throw(
            ArgumentError(
                "cannot reinterpret Tableau: backing bytes not divisible by sizeof($(U))",
            ),
        )
    end
    new_rows = (total_bytes) ÷ sizeof(U)
    if new_rows % 2 != 0
        throw(
            ArgumentError(
                "cannot reinterpret Tableau: resulting backing array length $(new_rows) is not even (cannot split into X and Z halves).",
            ),
        )
    end
    raw = vec(t.xzs)
    new_vec = Base.reinterpret(U, raw)
    new_xzs = reshape(new_vec, new_rows, cols)
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
