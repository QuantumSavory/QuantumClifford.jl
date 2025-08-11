
#=============================================================================#
# CAUTION: Keep in mind that the constants match the direction of the order.
# TODO: Make the parameters keyword arguments once support becomes available.
KA.@kernel inbounds = true unsafe_indices = true function kernel_mul!(
    mutable_phases::AbstractArray{<: Unsigned}, mutable_xzs::AbstractArray{T},
    @Const(const_xzs::AbstractArray{T}),
    multiplication_order::MultiplicationOrder,
    ::Val{phases}, ::Val{primary_axis}, ::Val{block_size}, ::Val{batch_size}
    ) where {
        T <: Unsigned, phases, primary_axis, block_size, batch_size
        }

    if primary_axis == primary_axis_rows
        j_mutable, begin_i = global_index(
            KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
            )
        stride_i = KA.@ndrange()[0x2]
    elseif primary_axis == primary_axis_qubits
        begin_i, j_mutable = global_index(
            KA.@index(Group, NTuple), KA.@groupsize(), KA.@index(Local, NTuple)
            )
        stride_i = KA.@ndrange()[0x1]
    end
    end_i = KA.@uniform (size(mutable_xzs, 0x1) >> 0x1)
    flag = KA.@uniform (size(const_xzs, 0x2) > 0x1)
    j_const = ifelse(flag, j_mutable, one(j_mutable))

    if phases
        low = zero(T)
        high = zero(T)
        phase_buffer = KA.@localmem DeviceUnsigned block_size
    end

    mutable_xzs = @view mutable_xzs[:, j_mutable]
    if multiplication_order == multiplication_order_left
        left = @view const_xzs[:, j_const]
        right = mutable_xzs
    elseif multiplication_order == multiplication_order_right
        left = mutable_xzs
        right = @view const_xzs[:, j_const]
    end

    for (i, _) in zip(begin_i : stride_i : end_i, one(batch_size) : batch_size)
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

        mutable_xzs[i] = x_new
        mutable_xzs[i + end_i] = z_new
    end

    if phases
        local_index = KA.@index(Local, Linear)
        phase_buffer[local_index] =
            ((count_ones(high) << 0x1) + count_ones(low)) & 0x3
        shared_memory_reduce!(
            reduce_sum!, local_index, Val(block_size), phase_buffer
            )

        if local_index == one(local_index)
            # CAUTION: This is sufficient since only atomicity is required.
            @atomic :monotonic mutable_phases[j_mutable] +=
                phase_buffer[local_index] & 0x3
            @atomic :monotonic mutable_phases[j_mutable] &= 0x3
        end
    end

end

# CAUTION: Requires either rows(const) == 1 or rows(const) == rows(mutable)
function device_mul!(
    mutable_phases::AbstractArray{<: Unsigned}, mutable_xzs::AbstractArray{T},
    const_phases::AbstractArray{<: Unsigned}, const_xzs::AbstractArray{T},
    multiplication_order::MultiplicationOrder;
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    )::Nothing where {
        T <: Unsigned, phase_B, primary_axis_E, block_SZ, batch_SZ
        }

    phase_B isa Bool && primary_axis_E isa PrimaryAxis &&
        block_SZ isa Integer && block_SZ > zero(block_SZ) &&
            batch_SZ isa Integer && batch_SZ > zero(batch_SZ) ||
                throw(ArgumentError(THROW_VALS))

    backend = KA.get_backend(mutable_xzs)

    if primary_axis_E == primary_axis_rows
        tile = (one(block_SZ), block_SZ)
        space = tessellate(
            (
                size(mutable_xzs, 0x2),
                cld(size(mutable_xzs, 0x1) >> 0x1, batch_SZ)
                ),
            tile
            )
    elseif primary_axis_E == primary_axis_qubits
        tile = (block_SZ, one(block_SZ))
        space = tessellate(
            (
                cld(size(mutable_xzs, 0x1) >> 0x1, batch_SZ),
                size(mutable_xzs, 0x2)
                ),
            tile
            )
    end

    if phase_B
        snippet! = kernel_snippet!(backend)
        @inbounds snippet!(
            snippet_mod_4_sum_phase!,
            mutable_phases, const_phases;
            ndrange = length(mutable_phases)
            )
    end
    mul! = kernel_mul!(backend)
    mul!(
        mutable_phases, mutable_xzs, const_xzs, multiplication_order,
        phases, primary_axis, block_size, batch_size;
        workgroupsize = tile, ndrange = space
        )

    return nothing

end
#=============================================================================#
