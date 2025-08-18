
#=============================================================================#
# TODO: Include the unsafe functions once the main package establishes them.
import QuantumClifford: canonicalize_rref!

# TODO: Implement support for (Mixed)Destabilizers.
function device_canonicalize_rref!(
    ph::AbstractArray{<: Unsigned}, xzs::AbstractArray{T},
    output_buffer::Union{Nothing, AbstractArray{<: Integer}} = nothing,
    bit_masks::Union{Nothing, AbstractArray{T}} = nothing;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    pauli_preferance::PauliPreferance = default_pauli_preferance,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )::Nothing where {T <: Unsigned}

    backend = KA.get_backend(xzs)

    if primary_axis == primary_axis_rows
        tile = (one(block_size), block_size)
        space = tessellate(
            (size(xzs, 2), cld(size(xzs, 1) >> 1, batch_size)),
            tile
            )
    elseif primary_axis == primary_axis_qubits
        tile = (block_size, one(block_size))
        space = tessellate(
            (cld(size(xzs, 1) >> 1, batch_size), size(xzs, 2)),
            tile
            )
    end

    # Utilised for loop management.
    length_xzs = size(xzs, 1) >> 1
    row_count = size(xzs, 2)
    toggle = false
    shrink_workspace = isnothing(bit_masks)

    # Required for safety whilst setting up for the proceeding iteration.
    mutex = create_mutex(backend)
    # Double buffered for present and preceeding/proceeding iteration.
    stride_fill = tracker_element_count
    tracker = similar(xzs, Csize_t, stride_fill << 1)
    fill!(tracker, typemax(Csize_t))
    # The pivot row tracker is initialised differently.
    begin_fill = Integer(tracker_content_swap_to)
    @inbounds fill!(
        (@view tracker[begin_fill : stride_fill : end]),
        row_count
        )

    if !isnothing(output_buffer)
        fill!(output_buffer, row_count)
    end

    bit_scan! = kernel_bit_scan!(backend)
    snippet! = kernel_snippet!(backend)
    mul_and_scan! = kernel_mul_and_scan!(backend)

    bit_scan!(
        tracker, mutex, xzs, bit_masks, primary_axis,
        Val(pauli_preferance), Val(sort_order_qubit_number),
        Val(block_size), Val(batch_size);
        workgroupsize = tile, ndrange = space
        )
    for _ in Base.OneTo(row_count)
        snippet!(
            snippet_swap_rows_prepare_tracker!,
            ph, xzs, tracker, toggle;
            ndrange = length_xzs
            )
        snippet!(
            snippet_set_row_phase_flag!,
            ph, xzs, tracker, toggle, pauli_preferance;
            ndrange = row_count
            )
        mul_and_scan!(
            ph, xzs, tracker, toggle, mutex, bit_masks, shrink_workspace,
            scan_side_lesser, multiplication_order, primary_axis,
            Val(pauli_preferance), Val(sort_order_qubit_number),
            Val(phases), Val(block_size), Val(batch_size);
            workgroupsize = tile, ndrange = space
            )
        # Switching the toggle is intentional.
        snippet!(
            snippet_track_pivot_canonicalize_rref!,
            output_buffer, tracker, !toggle;
            ndrange = 1
            )

        toggle = !toggle
    end

    # Remove all extraneous bits left behind by previous iterations.
    snippet!(snippet_mod_4_phase!, ph; ndrange = row_count)

    return nothing

end

# Tableau/AbstractStabilizer
@inline function canonicalize_rref!(
    state::DeviceUnionTableau,
    output_buffer::Union{Nothing, AbstractGPUArray{<: Integer}} = nothing,
    bit_masks::Union{Nothing, AbstractGPUArray{<: Unsigned}} = nothing;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    pauli_preferance::PauliPreferance = default_pauli_preferance,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    @boundscheck begin
        block_size > zero(block_size) && batch_size > zero(batch_size) ||
            throw(DomainError(THROW_PARAMETERS))
        if !isnothing(output_buffer)
            length(output_buffer) == 1 ||
                throw(ArgumentError(THROW_SIZE))
        end
        if !isnothing(bit_masks)
            length(bit_masks) == size(tab(state).xzs, 1) >> 1 ||
                throw(DimensionMismatch(THROW_SIZE))
        end
    end
    return do_canonicalize_rref!(
        state, output_buffer, bit_masks;
        multiplication_order = multiplication_order,
        pauli_preferance = pauli_preferance,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

end

@inline function do_canonicalize_rref!(
    state::DeviceUnionTableau,
    output_buffer::Union{Nothing, AbstractGPUArray{<: Integer}} = nothing,
    bit_masks::Union{Nothing, AbstractGPUArray{<: Unsigned}} = nothing;
    multiplication_order::MultiplicationOrder = default_multiplication_order,
    pauli_preferance::PauliPreferance = default_pauli_preferance,
    primary_axis::PrimaryAxis = default_primary_axis,
    phases::Bool = default_phases,
    block_size::Integer = default_block_size,
    batch_size::Integer = default_batch_size
    )

    if state isa AbstractStabilizer
        state_tab = tab(stabilizerview(state))
    elseif state isa Tableau
        state_tab = state
    end

    device_canonicalize_rref!(
        state_tab.phases, state_tab.xzs, output_buffer, bit_masks;
        multiplication_order = multiplication_order,
        pauli_preferance = pauli_preferance,
        primary_axis = primary_axis, phases = phases,
        block_size = block_size, batch_size = batch_size
        )

    if isnothing(output_buffer)
        return state
    else
        return state, output_buffer
    end

end
#=============================================================================#
