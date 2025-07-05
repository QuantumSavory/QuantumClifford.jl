# Single qubit Clifford group

# The group has a size of 24 in total
# As a stop gap solution, this is implemented using CliffordOperator.
# TODO: This could also be a place where interface is used, and we could use this stop gap implementation
# to test future implementation of the single qubit clifford group interface.
# TODO: Write test for this
struct SingleQubitCliffordGroup
    c::CliffordOperator
end

function Base.:(==)(l::SingleQubitCliffordGroup, r::SingleQubitCliffordGroup)
    l.c == r.c
end

Base.show(io::IO, s::SingleQubitCliffordGroup) = show(io, s.c)

Base.convert(::Type{SingleQubitCliffordGroup}, x::CliffordOperator) = SingleQubitCliffordGroup(x)

function Base.:(*)(l::SingleQubitCliffordGroup, r::SingleQubitCliffordGroup)
    SingleQubitCliffordGroup(l.c * r.c)
end

# This is on purpose to only include AbstractCliffordOperator
# so that the underlying CliffordOperator will always be within the single qubit clifford group
function Base.:(*)(l::Type{<:AbstractSingleQubitOperator}, r::SingleQubitCliffordGroup)
    SingleQubitCliffordGroup(l) * r
end

function Base.:(*)(l::SingleQubitCliffordGroup, r::Type{<:AbstractSingleQubitOperator})
    l * SingleQubitCliffordGroup(r)
end

function SingleQubitCliffordGroup(s::Type{O}) where O<:AbstractSingleQubitOperator
    return SingleQubitCliffordGroup(CliffordOperator(s))
end

Base.one(::Type{SingleQubitCliffordGroup}) = SingleQubitCliffordGroup(one(CliffordOperator, 1))

Base.inv(s::SingleQubitCliffordGroup) = SingleQubitCliffordGroup(inv(s.c))

apply!(s::AbstractStabilizer, self::SingleQubitCliffordGroup, indices; phases::Bool=true) = apply!(s, self.c, indices; phases)