# Works even when broadcasting on zero-dimensional arrays.
@inline u32(v) = map(x -> UInt32(x), v)

# Surprisingly, these do not already exist.
@inline get_pauli(t::Tableau, i::Integer) =
    PauliOperator((@view t.phases[i]), t.nqubits, (@view t.xzs[:, i]))
@inline get_pauli(s::AbstractStabilizer, i::Integer) =
    PauliOperator(
        (@view s.tab.phases[i]),
        s.tab.nqubits,
        (@view s.tab.xzs[:, i])
        )

@inline phases(t::Tableau) = t.phases
@inline phases(t::Tableau, i::Integer) = (@view t.phases[i])
@inline phases(s::AbstractStabilizer) = s.tab.phases
@inline phases(s::AbstractStabilizer, i::Integer) = (@view s.tab.phases[i])
@inline xzs(t::Tableau) = t.xzs
@inline xzs(t::Tableau, i::Integer) = (@view t.xzs[:, i])
@inline xzs(s::AbstractStabilizer) = s.tab.xzs
@inline xzs(s::AbstractStabilizer, i::Integer) = (@view s.tab.xzs[:, i])
