
#=============================================================================#
# CAUTION: Utilises unsafe_indices = true, hence demanding boundary validation.
KA.@kernel inbounds = true unsafe_indices = true function kernel_snippet!(
    f!::Function, arguments...
    )

    global_position = global_index(
        KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
        )
    f!(global_position, arguments...)

end

#==============================================================================
SNIPPET CATALOGUE
==============================================================================#

@inline function snippet_mod_4_sum_phase!(
    global_position::NTuple{N, Integer}, phases::AbstractArray{<: Unsigned},
    partner::Union{Unsigned, AbstractArray{<: Unsigned}}
    )::Nothing where {N}

    @inbounds begin
        i = global_position[1]
        if i <= length(phases)
            if partner isa Integer
                phases[i] = (phases[i] + partner) & 0x3
            elseif partner isa AbstractArray
                j = ifelse(length(partner) > 1, i, one(i))
                phases[i] = (phases[i] + partner[j]) & 0x3
            end
        end
    end
    return nothing

end

@inline function snippet_mod_4_phase!(
    global_position::NTuple{N, Integer}, phases::AbstractArray{<: Unsigned}
    )::Nothing where {N}

    @inbounds begin
        i = global_position[1]
        if i <= length(phases)
            phases[i] &= 0x3
        end
    end
    return nothing

end

@inline function snippet_track_pivot_canonicalize!(
    global_position::NTuple{N, Integer},
    output_buffer::Union{Nothing, AbstractArray{<: Integer}},
    tracker::AbstractArray{<: Unsigned}, toggle::Bool
    )::Nothing where {N}

    @inbounds begin
        current = KA.@uniform (
            ifelse(toggle, tracker_element_count, zero(Csize_t))
            )
        previous = KA.@uniform (
            ifelse(toggle, zero(Csize_t), tracker_element_count)
            )
        if global_position[1] == 1
            bit_type = tracker[current + Integer(tracker_content_bit_type)]
            row = tracker[previous + Integer(tracker_content_swap_to)]
            invalid = Integer(pauli_bit_invalid)

            if !isnothing(output_buffer)

                primary = Integer(pauli_bit_primary)
                secondary = Integer(pauli_bit_secondary)
                previous_bit_type =
                    tracker[previous + Integer(tracker_content_bit_type)]
                # Primary => Invalid
                if bit_type >= invalid && previous_bit_type == primary
                    output_buffer[1] = row
                    output_buffer[2] = row
                # Secondary => Invalid
                elseif bit_type >= invalid && previous_bit_type == secondary
                    output_buffer[2] = row
                # Invalid/Primary => Secondary
                elseif bit_type == secondary && previous_bit_type != secondary
                    output_buffer[1] = row
                end

            end

            row += ifelse(bit_type < invalid, one(row), zero(row))
            tracker[current + Integer(tracker_content_swap_to)] = row
        end
    end
    return nothing

end

@inline function snippet_track_pivot_canonicalize_rref!(
    global_position::NTuple{N, Integer},
    output_buffer::Union{Nothing, AbstractArray{<: Integer}},
    tracker::AbstractArray{<: Unsigned}, toggle::Bool
    )::Nothing where {N}

    @inbounds begin
        current = KA.@uniform (
            ifelse(toggle, tracker_element_count, zero(Csize_t))
            )
        previous = KA.@uniform (
            ifelse(toggle, zero(Csize_t), tracker_element_count)
            )
        if global_position[1] == 1
            bit_type = tracker[current + Integer(tracker_content_bit_type)]
            row = tracker[previous + Integer(tracker_content_swap_to)]
            invalid = Integer(pauli_bit_invalid)

            if !isnothing(output_buffer)
                previous_bit_type =
                    tracker[previous + Integer(tracker_content_bit_type)]
                # Valid => Invalid
                if bit_type >= invalid && previous_bit_type < invalid
                    output_buffer[1] = row - one(row)
                end
            end

            row -= ifelse(bit_type < invalid, one(row), zero(row))
            tracker[current + Integer(tracker_content_swap_to)] = row
        end
    end
    return nothing

end

@inline function snippet_swap_rows_prepare_tracker!(
    global_position::NTuple{N, Integer},
    phases::AbstractArray{<: Unsigned}, xzs::AbstractArray{<: Unsigned},
    tracker::AbstractArray{S}, toggle::Bool
    )::Nothing where {N, S <: Unsigned}

    @inbounds begin
        i = global_position[1]
        end_i = KA.@uniform (size(xzs, 1) >> 1)
        current = KA.@uniform (
            ifelse(toggle, tracker_element_count, zero(Csize_t))
            )
        next = KA.@uniform (
            ifelse(toggle, zero(Csize_t), tracker_element_count)
            )
        valid =
            tracker[current + Integer(tracker_content_bit_type)] <
                Integer(pauli_bit_invalid)
        row_A = tracker[current + Integer(tracker_content_swap_from)]
        row_B = tracker[current + Integer(tracker_content_swap_to)]

        if valid && row_A != row_B
            if i <= end_i
                temp_x = xzs[i, row_A]
                temp_z = xzs[i + end_i, row_A]
                xzs[i, row_A] = xzs[i, row_B]
                xzs[i + end_i, row_A] = xzs[i + end_i, row_B]
                xzs[i, row_B] = temp_x
                xzs[i + end_i, row_B] = temp_z
            end
            if i == one(i)
                temp_phase = phases[row_A]
                phases[row_A] = phases[row_B]
                phases[row_B] = temp_phase
            end
        end
        if i == one(i)
            tracker[next + Integer(tracker_content_index)] = typemax(S)
            tracker[next + Integer(tracker_content_bit_shift)] = typemax(S)
            tracker[next + Integer(tracker_content_bit_type)] = typemax(S)
            tracker[next + Integer(tracker_content_swap_from)] = typemax(S)
        end
    end
    return nothing

end

@inline function snippet_set_row_phase_flag!(
    global_position::NTuple{N, Integer},
    phases::AbstractArray{P}, xzs::AbstractArray{T},
    tracker::AbstractArray{<: Unsigned}, toggle::Bool,
    pauli_preference::PauliPreference
    )::Nothing where {N, P <: Unsigned, T <: Unsigned}

    @inbounds begin
        z_offset = KA.@uniform (size(xzs, 1) >> 1)
        end_rows = KA.@uniform (size(xzs, 2))
        current = KA.@uniform (
            ifelse(toggle, tracker_element_count, zero(Csize_t))
            )
        row = global_position[1]
        index = tracker[current + Integer(tracker_content_index)]
        bit_shift = tracker[current + Integer(tracker_content_bit_shift)]
        bit_type = tracker[current + Integer(tracker_content_bit_type)]

        if bit_type < Integer(pauli_bit_invalid) && row <= end_rows
            if bit_type == Integer(pauli_bit_primary)
                index += ifelse(
                    pauli_preference == pauli_preference_x,
                    zero(z_offset),
                    z_offset
                    )
                status = xzs[index, row] & bit_mask(bit_shift, T)
            elseif bit_type == Integer(pauli_bit_secondary)
                index += ifelse(
                    pauli_preference == pauli_preference_x,
                    z_offset,
                    zero(z_offset)
                    )
                status = xzs[index, row] & bit_mask(bit_shift, T)
            end
            phases[row] &= 0x3
            if status != zero(T)
                # Top bit is utilised as auxiliary storage.
                phases[row] |= top_bits(0x1, P)
            end
        end
    end
    return nothing

end
#=============================================================================#
