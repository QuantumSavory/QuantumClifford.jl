# Works even when broadcasting on zero-dimensional arrays.
@inline function u32(v)
    return map(x -> UInt32(x), v)
end

# Surprisingly, these do not already exist.
@inline function get_pauli(t::Tableau, i::Integer)
    return PauliOperator(
        (@view t.phases[i]),
        t.nqubits,
        (@view t.xzs[:, i])
        )
end
@inline function get_pauli(s::AbstractStabilizer, i::Integer)
    return PauliOperator(
        (@view s.tab.phases[i]),
        s.tab.nqubits,
        (@view s.tab.xzs[:, i])
        )
end

@inline function phases(t::Tableau)
    return t.phases
end
@inline function phases(t::Tableau, i::Integer)
    return (@view t.phases[i])
end
@inline function phases(s::AbstractStabilizer)
    return s.tab.phases
end
@inline function phases(s::AbstractStabilizer, i::Integer)
    return (@view s.tab.phases[i])
end

@inline function xzs(t::Tableau)
    return t.xzs
end
@inline function xzs(t::Tableau, i::Integer)
    return (@view t.xzs[:, i])
end
@inline function xzs(s::AbstractStabilizer)
    return s.tab.xzs
end
@inline function xzs(s::AbstractStabilizer, i::Integer)
    return (@view s.tab.xzs[:, i])
end
