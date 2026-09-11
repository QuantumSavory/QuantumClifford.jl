
#=============================================================================#
function test_KA_canonicalization(synchronize, AT, cache)::Nothing
    for (round, nqubits) in Iterators.product(
        Base.OneTo(min_rounds), qubit_counts
        )

    # Keep the memory usage sane.
    rows = min(nqubits, max_rows)
    axis = primary_axes[round]
    block = block_sizes[round]
    batch = batch_sizes[round]

    @cached cache begin

        # Stabilizer
        host_stabilizer = Stabilizer(random_tableau(max_rows, nqubits))
        device_stabilizer = adapt(AT, host_stabilizer)

        # Placeholders
        host_temp_stabilizer = zero(host_stabilizer)
        device_temp_stabilizer = adapt(AT, host_temp_stabilizer)

        # Important optional argument.
        if isodd(round)
            colindices = one(nqubits) : (one(nqubits) << 1) : nqubits
            xzs = host_stabilizer.tab.xzs
            bit_masks = AT(zeros(eltype(xzs), size(xzs, 1) >> 1))
            fill!(bit_masks, alternating_bit_mask(eltype(xzs)))
        else
            colindices = Base.OneTo(nqubits)
            bit_masks = nothing
        end

        i = rand(Base.OneTo(rows))

        # canonicalize
        @test begin
            copy_to!(host_temp_stabilizer, host_stabilizer)
            copy_to!(device_temp_stabilizer, device_stabilizer)

            host_output, host_pivot_x, host_pivot_z = canonicalize!(
                host_temp_stabilizer; ranks = true
                )
            device_buffer = AT(zeros(typeof(host_pivot_x), 2))
            device_output, device_buffer = canonicalize!(
                device_temp_stabilizer, device_buffer;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            result = equal_phases(host_output, device_output, i)
            result &= equal_xzs(host_output, device_output, i)
            result &= [host_pivot_x, host_pivot_z] == Array(device_buffer)
        end

        # canonicalize_rref
        @test begin
            copy_to!(host_temp_stabilizer, host_stabilizer)
            copy_to!(device_temp_stabilizer, device_stabilizer)

            host_output, host_pivot = canonicalize_rref!(
                host_temp_stabilizer, colindices
                )
            device_buffer = AT(zeros(typeof(host_pivot), 1))
            device_output, device_buffer = canonicalize_rref!(
                device_temp_stabilizer, device_buffer, bit_masks;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            result = equal_phases(host_output, device_output, i)
            result &= equal_xzs(host_output, device_output, i)
            result &= [host_pivot] == Array(device_buffer)
        end

        # canonicalize_gott
        # TODO: Implement.

    # Marks the end for @cached
    end

    # Marks the end for (round, nqubits)
    end

    return nothing
end
#=============================================================================#
