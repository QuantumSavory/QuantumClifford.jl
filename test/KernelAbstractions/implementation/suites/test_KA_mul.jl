
#=============================================================================#
# TODO: Revisit should the base package support more multiplication signatures.
# This must be done explicitly as they are not exported.
using QuantumClifford: mul_left!, mul_right!

function test_KA_mul(synchronize, AT, cache)::Nothing
    for (round, nqubits) in Iterators.product(
        Base.OneTo(max_rounds), qubit_counts
        )

    # Keep the memory usage sane.
    rows = min(nqubits, max_rows)
    axis = primary_axes[round]
    block = block_sizes[round]
    batch = batch_sizes[round]

    for mul! in (mul_left!, mul_right!)
    @cached cache begin

        # PauliOperator
        host_pauli_1 = random_pauli(nqubits)
        device_pauli_1 = adapt(AT, host_pauli_1)
        host_pauli_2 = random_pauli(nqubits)
        device_pauli_2 = adapt(AT, host_pauli_2)

        # Tableau
        host_tableau = random_tableau(rows, nqubits)
        device_tableau = adapt(AT, host_tableau)

        # Destabilizer
        host_destabilizer = Destabilizer(random_tableau(rows << 1, nqubits))
        device_destabilizer = adapt(AT, host_destabilizer)

        # Placeholders
        host_temp_pauli = zero(host_pauli_1)
        device_temp_pauli = adapt(AT, host_temp_pauli)
        host_temp_tableau = zero(host_tableau)
        device_temp_tableau = adapt(AT, host_temp_tableau)

        i = rand(Base.OneTo(rows))
        j = rand(Base.OneTo(rows))

        # Potential aliasing problem.
        @test begin
            copy_to!(host_temp_tableau, host_tableau)
            copy_to!(device_temp_tableau, device_tableau)

            host_output = mul!(host_temp_tableau, i, i)
            device_output = mul!(
                device_temp_tableau, i, i;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = equal_phases(host_output, device_output, i)
            flag &= equal_xzs(host_output, device_output, i)
        end

        # PauliOperator - PauliOperator
        @test begin
            copy_to!(host_temp_pauli, host_pauli_1)
            copy_to!(device_temp_pauli, device_pauli_1)

            host_output = mul!(host_temp_pauli, host_pauli_2)
            device_output = mul!(
                device_temp_pauli, device_pauli_2;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = host_output.phase == Array(device_output.phase)
            flag &= host_output.xz == Array(device_output.xz)
        end

        # PauliOperator - Tableau/AbstractStabilizer[i]
        @test begin
            copy_to!(host_temp_pauli, host_pauli_1)
            copy_to!(device_temp_pauli, device_pauli_1)

            host_output = mul!(host_temp_pauli, host_tableau, i)
            device_output = mul!(
                device_temp_pauli, device_tableau, i;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = host_output.phase == Array(device_output.phase)
            flag &= host_output.xz == Array(device_output.xz)
        end

        # Tableau/AbstractStabilizer - PauliOperator
        @test begin
            copy_to!(host_temp_tableau, host_tableau)
            copy_to!(device_temp_tableau, device_tableau)

            host_output = mul!(host_temp_tableau, host_pauli_1)
            device_output = mul!(
                device_temp_tableau, device_pauli_1;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = equal_phases(host_output, device_output, i)
            flag &= equal_xzs(host_output, device_output, i)
        end

        # Tableau/AbstractStabilizer - Tableau/AbstractStabilizer[i]
        @test begin
            copy_to!(host_temp_tableau, host_tableau)
            copy_to!(device_temp_tableau, device_tableau)

            host_output = mul!(host_temp_tableau, view_pauli(host_tableau, i))
            device_output = mul!(
                device_temp_tableau, device_tableau, i;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = equal_phases(host_output, device_output, i)
            flag &= equal_xzs(host_output, device_output, i)
        end

        # Tableau - Tableau/AbstractStabilizer
        @test begin
            copy_to!(device_temp_tableau, device_tableau)

            device_output = mul!(
                device_temp_tableau, device_tableau;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            # Rows become {+/-} x Identity since it is multiplied by itself.
            flag = reduce(|, view_phases(device_output)) & 0x1 == 0
            flag &= reduce(|, view_xzs(device_output, i)) == 0
        end

        # Tableau[i] - PauliOperator
        @test begin
            copy_to!(host_temp_pauli, view_pauli(host_tableau, i))
            copy_to!(device_temp_tableau, device_tableau)

            host_output = mul!(host_temp_pauli, host_pauli_1)
            device_output = mul!(
                device_temp_tableau, i, device_pauli_1;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = host_output.phase == Array(view_phases(device_output, i))
            flag &= host_output.xz == Array(view_xzs(device_output, i))
        end

        # Tableau[i] - Tableau/AbstractStabilizer[j]
        @test begin
            copy_to!(host_temp_tableau, host_tableau)
            copy_to!(device_temp_tableau, device_tableau)

            host_output = mul!(host_temp_tableau, i, host_tableau, j)
            device_output = mul!(
                device_temp_tableau, i, device_tableau, j;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = equal_phases(host_output, device_output, i)
            flag &= equal_xzs(host_output, device_output, i)
        end

        # Tableau/AbstractStabilizer[i] - Self[j]
        @test begin
            copy_to!(host_temp_tableau, host_tableau)
            copy_to!(device_temp_tableau, device_tableau)

            host_output = mul!(host_temp_tableau, i, j)
            device_output = mul!(
                device_temp_tableau, i, j;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = equal_phases(host_output, device_output, i)
            flag &= equal_xzs(host_output, device_output, i)
        end

        # (Mixed)Destabilizer[i] - Self[j]
        @test begin
            # There is no need to copy.
            host_output_1 = mul!(
                view_pauli(host_destabilizer, j),
                view_pauli(host_destabilizer, i);
                phases = Val(false)
                )
            offset = length(host_destabilizer.tab.phases) >> 1
            host_output_2 = mul!(
                view_pauli(host_destabilizer, i + offset),
                view_pauli(host_destabilizer, j + offset)
                )
            device_output = mul!(
                device_destabilizer, i, j;
                primary_axis = axis, block_size = block, batch_size = batch
                )
            synchronize()

            flag = host_output_1.phase == Array(view_phases(device_output, j))
            flag &= host_output_1.xz == Array(view_xzs(device_output, j))
            i += offset
            flag &= host_output_2.phase == Array(view_phases(device_output, i))
            flag &= host_output_2.xz == Array(view_xzs(device_output, i))
        end

    # Marks the end for @cached
    end
    # Marks the end for mul!
    end

    # Marks the end for (round, nqubits)
    end

    return nothing
end
#=============================================================================#
