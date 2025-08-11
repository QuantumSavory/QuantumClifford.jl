
#=============================================================================#
import QuantumClifford: canonicalize!

function device_canonicalize!(
    ph::AbstractArray{<: Unsigned}, xzs::AbstractArray{<: Unsigned},
    output_buffer::Union{Nothing, AbstractArray{S}} = nothing;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    )::Nothing where {
        S <: Integer, phase_B, primary_axis_E, block_SZ, batch_SZ
        }

    phase_B isa Bool && primary_axis_E isa PrimaryAxis &&
        block_SZ isa Integer && block_SZ > zero(block_SZ) &&
            batch_SZ isa Integer && batch_SZ > zero(batch_SZ) ||
                throw(ArgumentError(THROW_VALS))

    backend = KA.get_backend(xzs)

    if primary_axis_E == primary_axis_rows
        tile = (one(block_SZ), block_SZ)
        space = tessellate(
            (size(xzs, 0x2), cld(size(xzs, 0x1) >> 0x1, batch_SZ)),
            tile
            )
    elseif primary_axis_E == primary_axis_qubits
        tile = (block_SZ, one(block_SZ))
        space = tessellate(
            (cld(size(xzs, 0x1) >> 0x1, batch_SZ), size(xzs, 0x2)),
            tile
            )
    end

    # Utilised for loop management.
    length_xzs = size(xzs, 0x1) >> 0x1
    row_count = size(xzs, 0x2)
    toggle = false
    cycles_until_sync = default_scheduling_limit

    # Required for safety whilst setting up for the proceeding iteration.
    mutex = create_mutex(backend)
    # Double buffered for present and preceeding/proceeding iteration.
    stride_fill = tracker_element_count
    tracker = similar(xzs, Csize_t, stride_fill << 0x1)
    fill!(tracker, typemax(Csize_t))
    # The pivot row tracker is initialised differently.
    begin_fill = Integer(tracker_content_swap_to)
    @inbounds fill!(
        (@view tracker[begin_fill : stride_fill : end]),
        zero(Csize_t)
        )

    if !isnothing(output_buffer)
        fill!(output_buffer, zero(S))
    end

    bit_scan = kernel_bit_scan(backend)
    snippet! = kernel_snippet!(backend)
    mul_and_scan! = kernel_mul_and_scan!(backend)

    bit_scan(
        xzs, nothing, mutex, tracker,
        Val(sort_order_pauli_bit_prefer_x),
        primary_axis, block_size, batch_size;
        workgroupsize = tile, ndrange = space
        )
    snippet!(
        snippet_track_pivot_canonicalize!,
        output_buffer, tracker, toggle, sort_order_pauli_bit_prefer_x;
        ndrange = 0x1
        )
    for _ in one(row_count) : row_count
        snippet!(
            snippet_swap_rows_prepare_tracker!,
            ph, xzs, tracker, toggle;
            ndrange = length_xzs
            )
        snippet!(
            snippet_set_row_phase_flag!,
            ph, xzs, tracker, toggle;
            ndrange = row_count
            )
        mul_and_scan!(
            ph, xzs, multiplication_order, scan_side_greater,
            nothing, false, mutex, tracker, toggle,
            phases, Val(sort_order_pauli_bit_prefer_x),
            primary_axis, block_size, batch_size;
            workgroupsize = tile, ndrange = space
            )
        # Switching the toggle is intentional.
        snippet!(
            snippet_track_pivot_canonicalize!,
            output_buffer, tracker, !toggle, sort_order_pauli_bit_prefer_x;
            ndrange = 0x1
            )

        toggle = !toggle
        cycles_until_sync -= one(cycles_until_sync)
        if cycles_until_sync == zero(cycles_until_sync)
            KA.synchronize(backend)
            host_tracker = Array(tracker)
            offset = ifelse(toggle, tracker_element_count, 0x0)
            @inbounds continue_flag =
                host_tracker[offset + Integer(tracker_content_bit_type)] <
                    Integer(pauli_bit_invalid)
            if continue_flag
                cycles_until_sync = default_scheduling_limit
            else
                break
            end
        end
    end

    # Remove all extraneous bits left behind by previous iterations.
    snippet!(snippet_mod_4_phase!, ph; ndrange = row_count)

    return nothing

end

# Tableau
@inline function canonicalize!(
    tab::DeviceTableau,
    output_buffer::Union{Nothing, AbstractGPUArray{<: Integer}} = nothing;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        if !isnothing(output_buffer)
            length(output_buffer) == 0x2 ||
                throw(ArgumentError(THROW_SIZE))
        end
    end

    device_canonicalize!(
        tab.phases, tab.xzs, output_buffer;
        multiplication_order = multiplication_order,
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

    if isnothing(output_buffer)
        return tab
    else
        return tab, output_buffer
    end

end

# AbstractStabilizer
@inline function canonicalize!(
    state::DeviceAbstractStabilizer,
    output_buffer::Union{Nothing, AbstractGPUArray{<: Integer}} = nothing;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    phases::Val{phase_B} = Val(default_phases),
    primary_axis::Val{primary_axis_E} = Val(default_primary_axis),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, primary_axis_E, block_SZ, batch_SZ}

    @boundscheck begin
        if !isnothing(output_buffer)
            length(output_buffer) == 0x2 ||
                throw(ArgumentError(THROW_SIZE))
        end
    end

    if state isa Stabilizer
        upper = size(state.tab.xzs, 0x2)
        lower = one(upper)
    elseif state isa Destabilizer
        upper = size(state.tab.xzs, 0x2)
        lower = (upper >> one(upper)) + one(upper)
    elseif state isa MixedStabilizer
        upper = state.rank
        lower = one(upper)
    elseif state isa MixedDestabilizer
        upper = size(state.tab.xzs, 0x2)
        lower = (upper >> one(upper)) + one(upper)
        upper = (upper >> one(upper)) + state.rank
    end

    @inbounds device_canonicalize!(
        (@view state.tab.phases[lower : upper]),
        (@view state.tab.xzs[:, lower : upper]),
        output_buffer;
        multiplication_order = multiplication_order,
        phases = phases, primary_axis = primary_axis,
        block_size = block_size, batch_size = batch_size
        )

    if isnothing(output_buffer)
        return state
    else
        return state, output_buffer
    end

end
#=============================================================================#
