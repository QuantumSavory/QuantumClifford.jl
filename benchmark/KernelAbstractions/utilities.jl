# Works even when broadcasting on zero-dimensional arrays.
@inline u32(v) = map(x -> UInt32(x), v)

# By definition, the size of (unsigned) char is set to unity.
@inline bit_count(::Type{T}) where {T} = sizeof(T) * count_zeros(zero(Cuchar))
