
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

@inline function nqubits_pauli(size_MiB::Integer)
    # Each qubit requires 2 bits.
    return cld(size_MiB * MiB, 2)
end

@inline function nqubits_tableau(size_MiB::Integer)
    # Each qubit requires 2 bits.
    return round(Int, sqrt(cld(size_MiB * MiB, 2)), RoundUp)
end
#=============================================================================#
