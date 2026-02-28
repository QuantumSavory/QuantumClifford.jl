
#=============================================================================#
@inline function equal_phases(
    u::U, v::V
    ) where {
        U <: Union{Tableau, AbstractStabilizer},
        V <: Union{Tableau, AbstractStabilizer}
        }

    return Array(view_phases(u)) == Array(view_phases(v))

end
@inline function equal_phases(
    u::U, v::V, i::Integer
    ) where {
        U <: Union{Tableau, AbstractStabilizer},
        V <: Union{Tableau, AbstractStabilizer}
        }

    return Array(view_phases(u, i)) == Array(view_phases(v, i))

end

@inline function equal_xzs(
    u::U, v::V
    ) where {
        U <: Union{Tableau, AbstractStabilizer},
        V <: Union{Tableau, AbstractStabilizer}
        }

    return Array(view_xzs(u)) == Array(view_xzs(v))

end
@inline function equal_xzs(
    u::U, v::V, i::Integer
    ) where {
        U <: Union{Tableau, AbstractStabilizer},
        V <: Union{Tableau, AbstractStabilizer}
        }

    return Array(view_xzs(u, i)) == Array(view_xzs(v, i))

end
#=============================================================================#
