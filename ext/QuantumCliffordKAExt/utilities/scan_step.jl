
#=============================================================================#
# CAUTION: Output layout is (index, bit_shift, bit_type, break_flag).
@inline function scan_step(
    x_bits::T, z_bits::T, current_index::Integer,
    index::Integer, bit_shift::Integer, bit_type::PauliBit,
    ::Val{sort_order}
    ) where {T <: Unsigned, sort_order}

    current_shift_x = lowest_set_bit(x_bits)
    current_shift_z = lowest_set_bit(z_bits)

    if sort_order == sort_order_pauli_bit_prefer_x

        if current_shift_x < bit_count(T) && isless(
            (pauli_bit_x, current_index, current_shift_x),
            (bit_type, index, bit_shift)
            )

            output = (current_index, current_shift_x, pauli_bit_x, true)

        elseif current_shift_z < bit_count(T) && isless(
            (pauli_bit_z, current_index, current_shift_z),
            (bit_type, index, bit_shift)
            )

            output = (current_index, current_shift_z, pauli_bit_z, false)

        else

            output = (index, bit_shift, bit_type, false)

        end

    elseif sort_order == sort_order_pauli_bit_prefer_z

        # This reverses the ordering of the Pauli bits.
        if current_shift_z < bit_count(T) && isless(
            (xor(Integer(pauli_bit_z), 0x1), current_index, current_shift_z),
            (xor(Integer(bit_type), 0x1), index, bit_shift)
            )

            output = (current_index, current_shift_z, pauli_bit_z, true)

        elseif current_shift_x < bit_count(T) && isless(
            (xor(Integer(pauli_bit_x), 0x1), current_index, current_shift_x),
            (xor(Integer(bit_type), 0x1), index, bit_shift)
            )

            output = (current_index, current_shift_x, pauli_bit_x, false)

        else

            output = (index, bit_shift, bit_type, false)

        end

    elseif sort_order == sort_order_qubit_number_prefer_x

        candidate = min(
            (current_shift_x, pauli_bit_x),
            (current_shift_z, pauli_bit_z)
            )

        if @inbounds candidate[0x1] < bit_count(T) && isless(
            (current_index, candidate...),
            (index, bit_shift, bit_type)
            )

            output = (current_index, candidate..., true)

        else

            output = (index, bit_shift, bit_type, false)

        end

    elseif sort_order == sort_order_qubit_number_prefer_z

        # This reverses the ordering of the Pauli bits.
        candidate = min(
            (current_shift_x, xor(Integer(pauli_bit_x), 0x1)),
            (current_shift_z, xor(Integer(pauli_bit_z), 0x1))
            )

        @inbounds if candidate[0x1] < bit_count(T) && isless(
            (current_index, candidate...),
            (index, bit_shift, xor(Integer(bit_type), 0x1))
            )

            output = (current_index, candidate..., true)

        else

            output = (index, bit_shift, bit_type, false)

        end

    end

    return output
end
#=============================================================================#
