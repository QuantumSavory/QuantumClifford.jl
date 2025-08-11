# Works even when broadcasting on zero-dimensional arrays.
@inline function u32(v)
    return map(x -> UInt32(x), v)
end

# By definition, the size of (unsigned) char is set to unity.
@inline function bit_count(::Type{T}) where {T}
    return sizeof(T) * count_zeros(zero(Cuchar))
end
