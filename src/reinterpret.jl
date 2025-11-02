"""
Helpers to reinterpret packed bit representations between unsigned element types.

These methods provide convenient `Base.reinterpret(::Type{U}, x)` overloads for
`PauliOperator`, `Tableau`, `Stabilizer`, and `PauliFrame` so tests and code can
easily convert between `UInt8/UInt16/UInt32/...` packed layouts.
"""

import Base: reinterpret

"""
Pack a vector of Bool into a Vector{U} with little-endian bit order.

The bit with index 1 goes to the least-significant bit of the first word.
"""
function _pack_bits(bits::AbstractVector{Bool}, ::Type{U}) where {U<:Unsigned}
    n = length(bits)
    bits_per = sizeof(U) * 8
    nwords = (n + bits_per - 1) ÷ bits_per
    out = zeros(U, nwords)
    @inbounds for i in 1:n
        if bits[i]
            word = (i - 1) ÷ bits_per + 1
            pos = (i - 1) % bits_per
            out[word] |= (U(1) << pos)
        end
    end
    return out
end

"""
Reinterpret a single PauliOperator to use storage elements of type `U`.

This unpacks the X/Z bits and repacks them into words of type `U`.
"""
function reinterpret(::Type{U}, p::PauliOperator) where {U<:Unsigned}
    xs = xbit(p)
    zs = zbit(p)
    xs_packed = _pack_bits(xs, U)
    zs_packed = _pack_bits(zs, U)
    xz = vcat(xs_packed, zs_packed)
    return PauliOperator(p.phase[], p.nqubits, xz)
end

function reinterpret(::Type{U}, t::Tableau) where {U<:Unsigned}
    r, n = size(t)
    if n == 0 || r == 0
        return Tableau(copy(t.phases), n, zeros(U, 0, 0))
    end
    bits_per = sizeof(U) * 8
    words = (n + bits_per - 1) ÷ bits_per
    new_xzs = zeros(U, 2 * words, r)
    @inbounds for col in 1:r
        p = t[col]
        p2 = reinterpret(U, p)
        new_xzs[:, col] = p2.xz
    end
    return Tableau(copy(t.phases), n, new_xzs)
end

reinterpret(::Type{U}, s::Stabilizer) where {U<:Unsigned} = Stabilizer(reinterpret(U, tab(s)))

function reinterpret(::Type{U}, f::PauliFrame) where {U<:Unsigned}
    return PauliFrame(reinterpret(U, f.frame), copy(f.measurements))
end
