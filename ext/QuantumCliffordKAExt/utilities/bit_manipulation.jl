 
#=============================================================================#
# By definition, the size of (unsigned) char is set to unity.
@inline function bit_count(::Type{T}) where {T}
    return sizeof(T) * count_zeros(zero(Cuchar))
end

# CAUTION: Zero indexed shift, valid values are less than bit_count(T).
@inline function highest_set_bit(bit_field::T) where {T <: Unsigned}
    return T(leading_zeros(bit_field))
end

# CAUTION: Zero indexed shift, valid values are less than bit_count(T).
@inline function lowest_set_bit(bit_field::T) where {T <: Unsigned}
    return T(trailing_zeros(bit_field))
end

# CAUTION: Unsigned typing is intentional for branchless code generation.
@inline function top_bits(count::Unsigned, ::Type{T}) where {T <: Unsigned}
    return ~(~zero(T) >> count)
end

# CAUTION: Unsigned typing is intentional for branchless code generation.
@inline function bottom_bits(count::Unsigned, ::Type{T}) where {T <: Unsigned}
    return ~(~zero(T) << count)
end

# CAUTION: Unsigned typing is intentional for branchless code generation.
@inline function bit_mask(bit_shift::Unsigned, ::Type{T}) where {T <: Unsigned}
    return one(T) << bit_shift
end
#=============================================================================#
