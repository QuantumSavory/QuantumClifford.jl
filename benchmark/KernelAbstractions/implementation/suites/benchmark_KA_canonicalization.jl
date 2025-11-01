
#=============================================================================#
@inline function host_rref!(args...; kwargs...)::Nothing
    canonicalize_rref!(args...; kwargs...)
    return nothing
end

@inline function device_rref!(synchronize, args...; kwargs...)::Nothing
    canonicalize_rref!(args...; kwargs...)
    synchronize()
    return nothing
end

function benchmark_KA_canonicalize_rref!(
    synchronize, AT,
    index, parameters,
    host_time, extrapolate_host, host_fit,
    device_time, extrapolate_device, device_fit;
    enable_masks = false,
    phases = KAExt.default_phases,
    primary_axis = KAExt.default_primary_axis
    )

    error_flag = false

    # These definitions ought to be kept local to this scope.
    initial_fit = [0.0, 0.0, 0.0, 1.0]
    @inline function model(n, p)
        return @. p[1] + p[2] * n + p[3] * n^2 + p[4] * n^3
    end

    nqubits = nqubits_tableau(sizes_MiB[index])
    point_count = ifelse(
        include_threshold_point,
        index,
        index - one(index))
    fit_nqubits = nqubits_tableau.(sizes_MiB[Base.OneTo(point_count)])

    @cached cache begin
        # BenchmarkTools evaluates the setup block at the global scope.
        global host_stabilizer = Stabilizer(random_tableau(nqubits, nqubits))
        global device_stabilizer = adapt(AT, host_stabilizer)
        global host_temp = copy(host_stabilizer)
        global device_temp = copy(device_stabilizer)

        if enable_masks
            colindices = one(nqubits) : (one(nqubits) << 1) : nqubits
            xzs = host_stabilizer.tab.xzs
            bit_masks = AT(zeros(eltype(xzs), size(xzs, 1) >> 1))
            fill!(bit_masks, alternating_bit_mask(eltype(xzs)))
        else
            colindices = Base.OneTo(nqubits)
            bit_masks = nothing
        end

        if extrapolate_host[1]
            host_time[index] = model(nqubits, host_fit[1])
            time_span = Second(0)
        else
            # Trigger compilation before benchmarking.
            host_rref!(host_temp, colindices; phases = phases)
            start = now()
            host_time[index] = @belapsed host_rref!(
                $host_temp, $colindices; phases = $phases
                ) evals = evals samples = samples seconds = seconds setup = (
                    copy_to!(host_temp, host_stabilizer);
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
                device_rref!(
                    synchronize, device_temp, nothing, bit_masks;
                    phases = phases, primary_axis = primary_axis,
                    block_size = block, batch_size = batch
                    )
                start = now()
                device_time[index, parameters_index] = @belapsed device_rref!(
                        $synchronize, $device_temp, nothing, $bit_masks;
                        phases = $phases, primary_axis = $primary_axis,
                        block_size = $block, batch_size = $batch
                        ) evals = evals samples = samples seconds = seconds setup = (
                            @cached cache copy_to!(device_temp, device_stabilizer);
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

    host_stabilizer = nothing
    device_stabilizer = nothing
    host_temp = nothing
    device_temp = nothing
    unsafe_free!(cache)
    GC.gc(true)

    return error_flag

end

@inline function benchmark_KA_canonicalization(synchronize, AT, path)::Nothing
    for (phase, axis) in Iterators.product(phases, primary_axes)
        temp = path * "/canonicalization/phases_$phase" * "_$axis"
        run_benchmark(
            benchmark_KA_canonicalize_rref!,
            synchronize, AT,
            temp * "/canonicalize_rref",
            "Tableau canonicalization\n(phases_$phase, $axis)",
            "Tableau size (MiB)";
            enable_masks = false, phases = phase, primary_axis = axis
            )
        run_benchmark(
            benchmark_KA_canonicalize_rref!,
            synchronize, AT,
            temp * "/canonicalize_rref_masked",
            "Tableau canonicalization (masked)\n(phases_$phase, $axis)",
            "Tableau size (MiB)";
            enable_masks = true, phases = phase, primary_axis = axis
            )
    end
    return nothing
end
#=============================================================================#
