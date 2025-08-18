
#=============================================================================#
# CAUTION: Output layout is (index, bit_shift, bit_type, break_flag).
@inline function scan_step(
    x_bits::T, z_bits::T, current_index::Integer,
    index::Integer, bit_shift::Integer, bit_type::PauliBit,
    ::Val{pauli_preferance}, ::Val{sort_order}
    ) where {T <: Unsigned, pauli_preferance, sort_order}

    if pauli_preferance == pauli_preferance_x
        current_shift_primary = lowest_set_bit(x_bits)
        current_shift_secondary = lowest_set_bit(z_bits)
    elseif pauli_preferance == pauli_preferance_z
        current_shift_primary = lowest_set_bit(z_bits)
        current_shift_secondary = lowest_set_bit(x_bits)
    end

    if sort_order == sort_order_pauli_bit

        if current_shift_primary < bit_count(T) && isless(
            (pauli_bit_primary, current_index, current_shift_primary),
            (bit_type, index, bit_shift)
            )

            output = (
                current_index, current_shift_primary,
                pauli_bit_primary, true
                )

        elseif current_shift_secondary < bit_count(T) && isless(
            (pauli_bit_secondary, current_index, current_shift_secondary),
            (bit_type, index, bit_shift)
            )

            output = (
                current_index, current_shift_secondary,
                pauli_bit_secondary, false
                )

        else

            output = (index, bit_shift, bit_type, false)

        end

    elseif sort_order == sort_order_qubit_number

        candidate = min(
            (current_shift_primary, pauli_bit_primary),
            (current_shift_secondary, pauli_bit_secondary)
            )

        if @inbounds candidate[1] < bit_count(T) && isless(
            (current_index, candidate...),
            (index, bit_shift, bit_type)
            )

            output = (current_index, candidate..., true)

        else

            output = (index, bit_shift, bit_type, false)

        end

    end

    return output
end
#=============================================================================#
