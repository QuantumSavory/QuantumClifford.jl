# Works even when broadcasting on zero-dimensional arrays.
@inline u32(v) = map(x -> UInt32(x), v)
