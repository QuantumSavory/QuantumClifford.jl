# Works even when broadcasting on zero-dimensional arrays.
@inline u32(v) = map(x -> UInt32(x), v)

@inline bit_count(::Type{T}) where {T} = count_zeros(zero(T))
