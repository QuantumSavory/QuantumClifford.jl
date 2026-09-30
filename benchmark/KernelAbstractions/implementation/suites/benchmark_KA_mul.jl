
#=============================================================================#
# This must be done explicitly as they are not exported.
using QuantumClifford: mul_left!, mul_right!

@inline function host_mul!(args...; kwargs...)::Nothing
    mul_left!(args...; kwargs...)
    return nothing
end

@inline function device_mul!(synchronize, args...; kwargs...)::Nothing
    mul_left!(args...; kwargs...)
    synchronize()
    return nothing
end

function benchmark_KA_mul_pauli_pauli!(
    synchronize, AT,
    index, parameters,
    host_time, extrapolate_host, host_fit,
    device_time, extrapolate_device, device_fit;
    phases = KAExt.default_phases,
    primary_axis = KAExt.default_primary_axis
    )

    # These definitions ought to be kept local to this scope.
    initial_fit = [0.0, 1.0]
    @inline function model(n, p)
        return @. p[1] + p[2] * n
    end

    # Each qubit requires 2 bits.
    nqubits = nqubits_pauli(sizes_MiB[index])
    point_count = ifelse(
        include_threshold_point,
        index,
        index - one(index))
    fit_nqubits = nqubits_pauli.(sizes_MiB[Base.OneTo(point_count)])

    @cached cache begin
        # BenchmarkTools evaluates the setup block at the global scope.
        global host_pauli = random_pauli(nqubits)
        global device_pauli = adapt(AT, host_pauli)
        global host_temp = copy(host_pauli)
        global device_temp = copy(device_pauli)
        host_pauli_mul = random_pauli(nqubits)
        device_pauli_mul = adapt(AT, host_pauli_mul)

        if extrapolate_host[1]
            host_time[index] = model(nqubits, host_fit[1])
            time_span = Second(0)
        else
            # Trigger compilation before benchmarking.
            host_mul!(host_temp, host_pauli_mul; phases = Val(phases))
            start = now()
            host_time[index] = @belapsed host_mul!(
                $host_temp, $host_pauli_mul; phases = Val($phases)
                ) evals = evals samples = samples seconds = seconds setup = (
                    copy_to!(host_temp, host_pauli);
                    )
            time_span = now() - start
        end

        if time_span >= Second(extrapolation_threshold)
            if point_count >= length(initial_fit)
                fit_time = host_time[Base.OneTo(point_count)]
                fit = curve_fit(model, fit_nqubits, fit_time, initial_fit)

                extrapolate_host[1] = true
                host_fit[1] = fit.param

                if !include_threshold_point
                    host_time[index] = model(nqubits, host_fit[1])
                end
            elseif host_permit_simple_scaling
                fit = copy(initial_fit)
                fit[end] = host_time[index] / model(nqubits, fit)
                extrapolate_host[1] = true
                host_fit[1] = fit
            else
                return true
            end
        end

        for (parameters_index, (block, batch)) in pairs(parameters)
            if extrapolate_device[parameters_index]
                device_time[index, parameters_index] =
                    model(nqubits, device_fit[parameters_index])
                time_span = Second(0)
            else
                device_mul!(
                    synchronize, device_temp, device_pauli_mul;
                    phases = Val(phases), primary_axis = primary_axis,
                    block_size = block, batch_size = batch
                    )
                start = now()
                device_time[index, parameters_index] = @belapsed device_mul!(
                        $synchronize, $device_temp, $device_pauli_mul;
                        phases = Val($phases), primary_axis = $primary_axis,
                        block_size = $block, batch_size = $batch
                        ) evals = evals samples = samples seconds = seconds setup = (
                            @cached cache copy_to!(device_temp, device_pauli);
                            )
                time_span = now() - start
            end

            if time_span >= Second(extrapolation_threshold)
                if point_count >= length(initial_fit)
                    fit_time =
                        device_time[Base.OneTo(point_count), parameters_index]
                    fit = curve_fit(model, fit_nqubits, fit_time, initial_fit)

                    extrapolate_device[parameters_index] = true
                    device_fit[parameters_index] = fit.param

                    if !include_threshold_point
                        device_time[index, parameters_index] =
                            model(nqubits, device_fit[parameters_index])
                    end
                elseif device_permit_simple_scaling
                    fit = copy(initial_fit)
                    fit[end] = device_time[index, parameters_index] / model(nqubits, fit)
                    extrapolate_device[parameters_index] = true
                    device_fit[parameters_index] = fit
                else
                    return true
                end
            end
        end
    end

    host_pauli = nothing
    device_pauli = nothing
    host_temp = nothing
    device_temp = nothing
    unsafe_free!(cache)
    GC.gc(true)

    return false

end

function benchmark_KA_mul_tableau_pauli!(
    synchronize, AT,
    index, parameters,
    host_time, extrapolate_host, host_fit,
    device_time, extrapolate_device, device_fit;
    phases = KAExt.default_phases,
    primary_axis = KAExt.default_primary_axis
    )

    error_flag = false

    # These definitions ought to be kept local to this scope.
    initial_fit = [0.0, 0.0, 1.0]
    @inline function model(n, p)
        return @. p[1] + p[2] * n + p[3] * n^2
    end

    nqubits = nqubits_tableau(sizes_MiB[index])
    point_count = ifelse(
        include_threshold_point,
        index,
        index - one(index))
    fit_nqubits = nqubits_tableau.(sizes_MiB[Base.OneTo(point_count)])

    @cached cache begin
        # BenchmarkTools evaluates the setup block at the global scope.
        global host_tableau = random_tableau(nqubits, nqubits)
        global device_tableau = adapt(AT, host_tableau)
        global host_temp = copy(host_tableau)
        global device_temp = copy(device_tableau)
        host_pauli = random_pauli(nqubits)
        device_pauli = adapt(AT, host_pauli)

        if extrapolate_host[1]
            host_time[index] = model(nqubits, host_fit[1])
            time_span = Second(0)
        else
            # Trigger compilation before benchmarking.
            host_mul!(host_temp, host_pauli; phases = Val(phases))
            start = now()
            host_time[index] = @belapsed host_mul!(
                $host_temp, $host_pauli; phases = Val($phases)
                ) evals = evals samples = samples seconds = seconds setup = (
                    copy_to!(host_temp, host_tableau);
                    )
            time_span = now() - start
        end

        if time_span >= Second(extrapolation_threshold)
            if point_count >= length(initial_fit)
                fit_time = host_time[Base.OneTo(point_count)]
                fit = curve_fit(model, fit_nqubits, fit_time, initial_fit)

                extrapolate_host[1] = true
                host_fit[1] = fit.param

                if !include_threshold_point
                    host_time[index] = model(nqubits, host_fit[1])
                end
            elseif host_permit_simple_scaling
                fit = copy(initial_fit)
                fit[end] = host_time[index] / model(nqubits, fit)
                extrapolate_host[1] = true
                host_fit[1] = fit
            else
                return true
            end
        end

        for (parameters_index, (block, batch)) in pairs(parameters)
            if extrapolate_device[parameters_index]
                device_time[index, parameters_index] =
                    model(nqubits, device_fit[parameters_index])
                time_span = Second(0)
            else
                device_mul!(
                    synchronize, device_temp, device_pauli;
                    phases = Val(phases), primary_axis = primary_axis,
                    block_size = block, batch_size = batch
                    )
                start = now()
                device_time[index, parameters_index] = @belapsed device_mul!(
                        $synchronize, $device_temp, $device_pauli;
                        phases = Val($phases), primary_axis = $primary_axis,
                        block_size = $block, batch_size = $batch
                        ) evals = evals samples = samples seconds = seconds setup = (
                            @cached cache copy_to!(device_temp, device_tableau);
                            )
                time_span = now() - start
            end

            if time_span >= Second(extrapolation_threshold)
                if point_count >= length(initial_fit)
                    fit_time =
                        device_time[Base.OneTo(point_count), parameters_index]
                    fit = curve_fit(model, fit_nqubits, fit_time, initial_fit)

                    extrapolate_device[parameters_index] = true
                    device_fit[parameters_index] = fit.param

                    if !include_threshold_point
                        device_time[index, parameters_index] =
                            model(nqubits, device_fit[parameters_index])
                    end
                elseif device_permit_simple_scaling
                    fit = copy(initial_fit)
                    fit[end] = device_time[index, parameters_index] / model(nqubits, fit)
                    extrapolate_device[parameters_index] = true
                    device_fit[parameters_index] = fit
                else
                    return true
                end
            end
        end
    end

    host_tableau = nothing
    device_tableau = nothing
    host_temp = nothing
    device_temp = nothing
    unsafe_free!(cache)
    GC.gc(true)

    return error_flag

end

@inline function benchmark_KA_mul(synchronize, AT, path)::Nothing
    for (phase, axis) in Iterators.product(phases, primary_axes)
        temp = path * "/mul/phases_$phase" * "_$axis"
        run_benchmark(
            benchmark_KA_mul_pauli_pauli!,
            synchronize, AT,
            temp * "/pauli_pauli",
            "Pauli-Pauli multiplication\n(phases_$phase, $axis)",
            "Pauli operator size (MiB)";
            phases = phase, primary_axis = axis
            )
        run_benchmark(
            benchmark_KA_mul_tableau_pauli!,
            synchronize, AT,
            temp * "/tableau_pauli",
            "Tableau-Pauli multiplication\n(phases_$phase, $axis)",
            "Tableau size (MiB)";
            phases = phase, primary_axis = axis
            )
    end
    return nothing
end
#=============================================================================#
