# This must be done explicitly as they are not exported.
using QuantumClifford: mul_left!, mul_right!

@inline host_f(x, y; phases::Val{phase_B} = Val(true)) where {phase_B} =
    mul_left!(x, y; phases = phases)

@inline function device_f(
    x, y, synchronize;
    phases::Val{phase_B} = Val(true),
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {phase_B, block_SZ, batch_SZ}

    mul_left!(
        x, y;
        phases = phases, block_size = block_size, batch_size = batch_size
        )
    synchronize()

end

@inline function benchmark_KA_mul_leftright(
    AT, synchronize, path;
    phases::Val{phase_B} = Val(true)
    ) where {phase_B}

    host_time = zeros(Float64, length(n_MiB))
    device_time = zeros(Float64, length(batch_sizes), length(n_MiB))

    # Keep the memory usage sane.
    cache = AllocCache()
    for (i, n) in enumerate(n_MiB)
        @cached cache begin
            # Each qubit requires 2 bits.
            h_p = PauliOperator(
                zeros(Cuchar), n * MiB >> 1,
                zeros(UInt, cld(n * MiB, bit_count(UInt)))
                )
            d_p = PauliOperator(AT(u32(h_p.phase)), h_p.nqubits, AT(h_p.xz))
            synchronize()
            # Trigger compilation before benchmarking.
            host_f(h_p, h_p; phases = phases)
            host_time[i] = @belapsed host_f(
                $h_p, $h_p; phases = $phases
                ) evals = evals samples = samples seconds = seconds
            for (j, size) in enumerate(batch_sizes)
                device_f(
                    d_p, d_p, synchronize;
                    phases = phases, batch_size = Val(size)
                    )
                device_time[j, i] =
                    @belapsed device_f(
                        $d_p, $d_p, $synchronize;
                        phases = $phases, batch_size = Val($size)
                        ) evals = evals samples = samples seconds = seconds
            end
        end
    end
    unsafe_free!(cache)

    device_cat = [device_time[i, :] for i = 1 : length(batch_sizes)]
    title =
    "Performance uplift - multiplication (phases = $phase_B)\n\
    Host threads = " * string(Threads.nthreads()) * " / " *
    string(Sys.CPU_THREADS) * ", Device block size = $default_block_size"
    xlabel = "Pauli operator size (MiB)"
    label = hcat(("Device - batch size = " .* string.(batch_sizes))..., "Host")

    path = "$path/mul_leftright"
    mkpath(path)

    plot(
        n_MiB, 10^3 .* hcat(device_cat..., host_time);
        shape = :circle, xticks = n_MiB, xscale = :log2, yscale = :log10,
        title = title, label = label, xlabel = xlabel, ylabel = "Runtime (ms)",
        background_color = :transparent
        )
    savefig("$path/runtime_phase_$phase_B.$format")

    plot(
        n_MiB, map(x -> host_time ./ x, device_cat);
        shape = :circle, xticks = n_MiB, xscale = :log2, title = title,
        label = hcat(label[1 : end - 1]...), xlabel = xlabel,
        ylabel = "Ratio (host/device)", background_color = :transparent
        )
    savefig("$path/ratio_phase_$phase_B.$format")

end
