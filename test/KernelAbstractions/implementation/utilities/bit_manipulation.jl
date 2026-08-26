
#=============================================================================#
# Repeats 0x55 to fill out all the bits in the given type.
@inline function alternating_bit_mask(::Type{T}) where {T <: Unsigned}
    counter = count_zeros(zero(T)) >> one(T)
    pattern = one(T)
    while counter > one(counter)
        pattern |= pattern << counter
        counter >>= one(counter)
    end
    return pattern
end
#=============================================================================#
