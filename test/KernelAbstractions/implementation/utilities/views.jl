
#=============================================================================#
# Avoid scalar indexing to combat complaints from the device backend(s).
@inline function view_pauli(t::Union{Tableau, AbstractStabilizer}, i::Integer)
    t_tab = tab(t)
    return PauliOperator(
        (@view t_tab.phases[i]),
        t_tab.nqubits,
        (@view t_tab.xzs[:, i])
        )
end

@inline function view_phases(t::Union{Tableau, AbstractStabilizer})
    return tab(t).phases
end
@inline function view_phases(t::Union{Tableau, AbstractStabilizer}, i::Integer)
    return (@view tab(t).phases[i])
end

@inline function view_xzs(t::Union{Tableau, AbstractStabilizer})
    return tab(t).xzs
end
@inline function view_xzs(t::Union{Tableau, AbstractStabilizer}, i::Integer)
    return (@view tab(t).xzs[:, i])
end
#=============================================================================#
