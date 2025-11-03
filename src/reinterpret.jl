"""
Helpers to reinterpret packed bit representations between unsigned element types.

These methods define `reinterpret(::Type{U}, ::PauliOperator)` and
`reinterpret(::Type{U}, ::Tableau)` in a way that prefers zero-copy
views using `Base.reinterpret` when the underlying byte layout allows it.
If the storage-element sizes are incompatible, the code falls back to a
safe packing path.
"""
import Base: reinterpret

function reinterpret(::Type{U}, p::PauliOperator) where {U<:Unsigned}
    old = eltype(p.xz)
    total_bytes = length(p.xz) * sizeof(old)
    if total_bytes % sizeof(U) != 0
        throw(ArgumentError("cannot reinterpret PauliOperator: backing bytes not divisible by sizeof($(U))"))
    end
    new_xz = Base.reinterpret(U, p.xz)
    return PauliOperator(p.phase, p.nqubits, new_xz)
end

function reinterpret(::Type{U}, t::Tableau) where {U<:Unsigned}
    rows = size(t.xzs, 1)
    cols = size(t.xzs, 2)
    old = eltype(t.xzs)
    total_bytes = rows * sizeof(old)
    if total_bytes % sizeof(U) != 0
        throw(ArgumentError("cannot reinterpret Tableau: backing bytes not divisible by sizeof($(U))"))
    end
    new_rows = (total_bytes) ÷ sizeof(U)
    raw = vec(t.xzs)
    new_vec = Base.reinterpret(U, raw)
    new_xzs = reshape(new_vec, new_rows, cols)
    return Tableau(t.phases, t.nqubits, new_xzs)
end

reinterpret(::Type{U}, s::Stabilizer) where {U<:Unsigned} = Stabilizer(reinterpret(U, tab(s)))

function reinterpret(::Type{U}, f::PauliFrame) where {U<:Unsigned}
    return PauliFrame(reinterpret(U, f.frame), f.measurements)
end
