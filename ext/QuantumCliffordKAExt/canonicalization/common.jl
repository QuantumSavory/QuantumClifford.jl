
#=============================================================================#
# TODO: Make the parameters keyword arguments once support becomes available.
KA.@kernel inbounds = true unsafe_indices = true function kernel_bit_scan!(
    tracker::AbstractArray{S}, mutex::AbstractMutex,
    @Const(xzs::AbstractArray{T}),
    @Const(bit_masks::Union{Nothing, AbstractArray{T}}),
    @Const(primary_axis::PrimaryAxis),
    ::Val{pauli_preference}, ::Val{sort_order},
    ::Val{block_size}, ::Val{batch_size}
    ) where {
        S <: Unsigned, T <: Unsigned,
        pauli_preference, sort_order, block_size, batch_size
        }

    if primary_axis == primary_axis_rows
        j, begin_i = global_index(
            KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
            )
        stride_i = KA.@ndrange()[2]
    elseif primary_axis == primary_axis_qubits
        begin_i, j = global_index(
            KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
            )
        stride_i = KA.@ndrange()[1]
    end
    end_i = KA.@uniform (size(xzs, 1) >> 1)

    index::S = typemax(S)
    bit_shift::S = typemax(S)
    bit_type = pauli_bit_invalid
    index_buffer = KA.@localmem S block_size
    bit_shift_buffer = KA.@localmem S block_size
    bit_type_buffer = KA.@localmem DeviceUnsigned block_size

    scan_target = @view xzs[:, j]

    for (i, _) in zip(begin_i : stride_i : end_i, Base.OneTo(batch_size))
        x_bits = scan_target[i]
        z_bits = scan_target[i + end_i]
        if !isnothing(bit_masks)
            mask = bit_masks[i]
            x_bits &= mask
            z_bits &= mask
        end

        index, bit_shift, bit_type, break_flag = scan_step(
            x_bits, z_bits, i, index, bit_shift, bit_type,
            Val(pauli_preference), Val(sort_order)
            )
        break_flag && break
    end

    local_index = DeviceUnsigned(KA.@index(Local, Linear))
    bit_shift_buffer[local_index] = bit_shift
    index_buffer[local_index] = index
    bit_type_buffer[local_index] = Integer(bit_type)

    if sort_order == sort_order_pauli_bit

        shared_memory_reduce!(
            reduce_lexicographic_min!,
            local_index, Val(DeviceUnsigned(block_size)),
            bit_type_buffer, index_buffer, bit_shift_buffer
            )

    elseif sort_order == sort_order_qubit_number

        shared_memory_reduce!(
            reduce_lexicographic_min!,
            local_index, Val(DeviceUnsigned(block_size)),
            index_buffer, bit_shift_buffer, bit_type_buffer
            )

    end

    if local_index == one(local_index)

        enumless_bit_type = bit_type_buffer[local_index]
        # Avoid mutex contention if there is no valid contribution.
        if enumless_bit_type < Integer(pauli_bit_invalid)
            index = index_buffer[local_index]
            bit_shift = bit_shift_buffer[local_index]

            lock_mutex!(
                mutex,
                (@view tracker[Integer(tracker_content_index)]),
                (@view tracker[Integer(tracker_content_bit_shift)]),
                (@view tracker[Integer(tracker_content_bit_type)]),
                (@view tracker[Integer(tracker_content_swap_from)])
                )

            temp_index = tracker[Integer(tracker_content_index)]
            temp_bit_shift = tracker[Integer(tracker_content_bit_shift)]
            temp_bit_type = tracker[Integer(tracker_content_bit_type)]
            temp_row = tracker[Integer(tracker_content_swap_from)]

            if sort_order == sort_order_pauli_bit

                lower_than_current = isless(
                    (enumless_bit_type, index, bit_shift, j),
                    (temp_bit_type, temp_index, temp_bit_shift, temp_row)
                    )

            elseif sort_order == sort_order_qubit_number

                lower_than_current = isless(
                    (index, bit_shift, enumless_bit_type, j),
                    (temp_index, temp_bit_shift, temp_bit_type, temp_row)
                    )

            end

            if lower_than_current
                tracker[Integer(tracker_content_index)] = index
                tracker[Integer(tracker_content_bit_shift)] = bit_shift
                tracker[Integer(tracker_content_bit_type)] = enumless_bit_type
                tracker[Integer(tracker_content_swap_from)] = j
            end

            unlock_mutex!(
                mutex,
                (@view tracker[Integer(tracker_content_index)]),
                (@view tracker[Integer(tracker_content_bit_shift)]),
                (@view tracker[Integer(tracker_content_bit_type)]),
                (@view tracker[Integer(tracker_content_swap_from)])
                )
        end

    end

end

# CAUTION: Keep in mind that the constants match the direction of the order.
# TODO: Make the parameters keyword arguments once support becomes available.
KA.@kernel inbounds = true unsafe_indices = true function kernel_mul_and_scan!(
    ph::AbstractArray{P}, xzs::AbstractArray{T},
    tracker::AbstractArray{S}, toggle::Bool, mutex::AbstractMutex,
    @Const(bit_masks::Union{Nothing, AbstractArray{T}}),
    @Const(shrink_workspace::Bool), @Const(scan_side::ScanSide),
    @Const(multiplication_order::MultiplicationOrder),
    @Const(primary_axis::PrimaryAxis),
    ::Val{pauli_preference}, ::Val{sort_order},
    ::Val{phases}, ::Val{block_size}, ::Val{batch_size}
    ) where {
        S <: Unsigned, P <: Unsigned, T <: Unsigned,
        pauli_preference, sort_order, phases, block_size, batch_size
        }

    if primary_axis == primary_axis_rows
        j, begin_i = global_index(
            KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
            )
        stride_i = KA.@ndrange()[2]
    elseif primary_axis == primary_axis_qubits
        begin_i, j = global_index(
            KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
            )
        stride_i = KA.@ndrange()[1]
    end
    end_i = KA.@uniform (size(xzs, 1) >> 1)

    local_index = DeviceUnsigned(KA.@index(Local, Linear))
    # Layout is [index, row].
    shared_parameters = KA.@localmem S 2
    # Dense storage as flags in a bit field.
    shared_flags = KA.@localmem DeviceUnsigned 1
    if local_index == one(local_index)

        shared_flags[1] = ifelse(
            begin_i == one(begin_i),
            Integer(bit_field_flag_leader),
            zero(DeviceUnsigned)
            )

        current = ifelse(toggle, tracker_element_count, zero(Csize_t))
        temp_index = tracker[current + Integer(tracker_content_index)]
        temp_bit_type = tracker[current + Integer(tracker_content_bit_type)]
        temp_row = tracker[current + Integer(tracker_content_swap_to)]

        new_begin_i = ifelse(
            shrink_workspace,
            begin_i + temp_index - one(S),
            begin_i + zero(S)
            )

        # Opt for an early exit if any of these conditions fail.
        continue_flag =
            temp_row != j && begin_i <= new_begin_i <= end_i &&
                temp_bit_type < Integer(pauli_bit_invalid)

        if continue_flag
            shared_parameters[1] = temp_index
            shared_parameters[2] = temp_row

            if scan_side == scan_side_lesser
                flags = ifelse(
                    j < temp_row,
                    Integer(bit_field_flag_scan),
                    zero(DeviceUnsigned)
                    )
            elseif scan_side == scan_side_greater
                flags = ifelse(
                    j > temp_row,
                    Integer(bit_field_flag_scan),
                    zero(DeviceUnsigned)
                    )
            end

            # Top bit holds the multiplication flag.
            flags |= ifelse(
                ph[j] & top_bits(0x1, P) != zero(P),
                Integer(bit_field_flag_multiply),
                zero(DeviceUnsigned)
                )

            shared_flags[1] |= flags
        end

    end
    KA.@synchronize()
    flags = shared_flags[1]

    if flags & (
        Integer(bit_field_flag_scan) | Integer(bit_field_flag_multiply)
        ) != zero(DeviceUnsigned)

    if phases
        low = zero(T)
        high = zero(T)
        phase_buffer = KA.@localmem DeviceUnsigned block_size
    end

    index::S = typemax(S)
    bit_shift::S = typemax(S)
    bit_type = pauli_bit_invalid
    index_buffer = KA.@localmem S block_size
    bit_shift_buffer = KA.@localmem S block_size
    bit_type_buffer = KA.@localmem DeviceUnsigned block_size

    new_begin_i = ifelse(
        shrink_workspace,
        begin_i + shared_parameters[1] - one(S),
        begin_i + zero(S)
        )

    if begin_i <= new_begin_i <= end_i

    begin_i = new_begin_i

    # Equivalent to bit_scan.
    if flags & (
        Integer(bit_field_flag_scan) | Integer(bit_field_flag_multiply)
        ) == Integer(bit_field_flag_scan)

    scan_target = @view xzs[:, j]

    for (i, _) in zip(begin_i : stride_i : end_i, Base.OneTo(batch_size))
        x_bits = scan_target[i]
        z_bits = scan_target[i + end_i]
        if !isnothing(bit_masks)
            mask = bit_masks[i]
            x_bits &= mask
            z_bits &= mask
        end

        index, bit_shift, bit_type, break_flag = scan_step(
            x_bits, z_bits, i, index, bit_shift, bit_type,
            Val(pauli_preference), Val(sort_order)
            )
        break_flag && break
    end

    else

    write_xzs = @view xzs[:, j]
    if multiplication_order == multiplication_order_left
        left = @view xzs[:, shared_parameters[2]]
        right = write_xzs
    elseif multiplication_order == multiplication_order_right
        left = write_xzs
        right = @view xzs[:, shared_parameters[2]]
    end

    # Equivalent to mul!.
    if flags & (
        Integer(bit_field_flag_scan) | Integer(bit_field_flag_multiply)
        ) == Integer(bit_field_flag_multiply)

    for (i, _) in zip(begin_i : stride_i : end_i, Base.OneTo(batch_size))
        x_left = left[i]
        z_left = left[i + end_i]
        x_right = right[i]
        z_right = right[i + end_i]

        x_new = xor(x_left, x_right)
        z_new = xor(z_left, z_right)
        if phases
            xl_zr = x_left & z_right
            merged = xor(xl_zr, z_left & x_right)
            high = xor(high, xor(low, x_new, z_new, xl_zr) & merged)
            low = xor(low, merged)
        end

        write_xzs[i] = x_new
        write_xzs[i + end_i] = z_new
    end

    # Equivalent to joint mul! and a modified bit_scan.
    else

    for (i, _) in zip(begin_i : stride_i : end_i, Base.OneTo(batch_size))
        x_left = left[i]
        z_left = left[i + end_i]
        x_right = right[i]
        z_right = right[i + end_i]

        x_new = xor(x_left, x_right)
        z_new = xor(z_left, z_right)
        if phases
            xl_zr = x_left & z_right
            merged = xor(xl_zr, z_left & x_right)
            high = xor(high, xor(low, x_new, z_new, xl_zr) & merged)
            low = xor(low, merged)
        end

        write_xzs[i] = x_new
        write_xzs[i + end_i] = z_new

        if !isnothing(bit_masks)
            mask = bit_masks[i]
            x_new &= mask
            z_new &= mask
        end

        index, bit_shift, bit_type = scan_step(
            x_new, z_new, i, index, bit_shift, bit_type,
            Val(pauli_preference), Val(sort_order)
            )

    end

    # Marks the end for flags == Integer(bit_field_flag_multiply)/else
    end

    # Marks the end for flags == Integer(bit_field_flag_scan)/else
    end

    # Marks the end for begin_i <= new_begin_i <= end_i
    end

    # Multiplication phase reduction and update.
    if phases
        if flags & Integer(bit_field_flag_multiply) != zero(DeviceUnsigned)
            phase_buffer[local_index] =
                ((count_ones(high) << 1) + count_ones(low)) & 0x3
            shared_memory_reduce!(
                reduce_sum!,
                local_index, Val(DeviceUnsigned(block_size)),
                phase_buffer
                )

            if local_index == one(local_index)
                temp_ph = phase_buffer[local_index]
                if flags & Integer(bit_field_flag_leader) !=
                    zero(DeviceUnsigned)
                    temp_ph += ph[shared_parameters[2]]
                end
                # Avoid expensive operations when they are not required.
                if stride_i > block_size
                    # CAUTION: This memory order is sufficient.
                    @atomic :monotonic ph[j] += temp_ph & 0x3
                    # CAUTION: Avoid nullifying the multiplication flag bit!
                    @atomic :monotonic ph[j] &= 0x3 | top_bits(0x1, P)
                else
                    ph[j] += temp_ph & 0x3
                    ph[j] &= 0x3
                end
            end
        end
    end

    # Pauli bit scan reduction and update.
    if flags & Integer(bit_field_flag_scan) != zero(DeviceUnsigned)
        bit_shift_buffer[local_index] = bit_shift
        index_buffer[local_index] = index
        bit_type_buffer[local_index] = Integer(bit_type)

        if sort_order == sort_order_pauli_bit

            shared_memory_reduce!(
                reduce_lexicographic_min!,
                local_index, Val(DeviceUnsigned(block_size)),
                bit_type_buffer, index_buffer, bit_shift_buffer
                )

        elseif sort_order == sort_order_qubit_number

            shared_memory_reduce!(
                reduce_lexicographic_min!,
                local_index, Val(DeviceUnsigned(block_size)),
                index_buffer, bit_shift_buffer, bit_type_buffer
                )

        end

        if local_index == one(local_index)

            enumless_bit_type = bit_type_buffer[local_index]
            # Avoid mutex contention if there is no valid contribution.
            if enumless_bit_type < Integer(pauli_bit_invalid)
                index = index_buffer[local_index]
                bit_shift = bit_shift_buffer[local_index]
                next = ifelse(toggle, zero(Csize_t), tracker_element_count)

                lock_mutex!(
                    mutex,
                    (@view tracker[next + Integer(tracker_content_index)]),
                    (@view tracker[next + Integer(tracker_content_bit_shift)]),
                    (@view tracker[next + Integer(tracker_content_bit_type)]),
                    (@view tracker[next + Integer(tracker_content_swap_from)])
                    )

                temp_index = tracker[next + Integer(tracker_content_index)]
                temp_bit_shift =
                    tracker[next + Integer(tracker_content_bit_shift)]
                temp_bit_type =
                    tracker[next + Integer(tracker_content_bit_type)]
                temp_row = tracker[next + Integer(tracker_content_swap_from)]

                if sort_order == sort_order_pauli_bit

                    lower_than_current = isless(
                        (enumless_bit_type, index, bit_shift, j),
                        (temp_bit_type, temp_index, temp_bit_shift, temp_row)
                        )

                elseif sort_order == sort_order_qubit_number

                    lower_than_current = isless(
                        (index, bit_shift, enumless_bit_type, j),
                        (temp_index, temp_bit_shift, temp_bit_type, temp_row)
                        )

                end

                if lower_than_current
                    tracker[next + Integer(tracker_content_index)] = index
                    tracker[next + Integer(tracker_content_bit_shift)] =
                        bit_shift
                    tracker[next + Integer(tracker_content_bit_type)] =
                        enumless_bit_type
                    tracker[next + Integer(tracker_content_swap_from)] = j
                end


                unlock_mutex!(
                    mutex,
                    (@view tracker[next + Integer(tracker_content_index)]),
                    (@view tracker[next + Integer(tracker_content_bit_shift)]),
                    (@view tracker[next + Integer(tracker_content_bit_type)]),
                    (@view tracker[next + Integer(tracker_content_swap_from)])
                    )
            end

        end
    end

    # Marks the end for flags & (...) != zero(DeviceUnsigned)
    end

end
#=============================================================================#
