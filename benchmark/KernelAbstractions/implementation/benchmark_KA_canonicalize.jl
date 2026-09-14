# This must be done explicitly as they are not exported.
using QuantumClifford: canonicalize!, canonicalize_rref!, Stabilizer, MixedDestabilizer, random_stabilizer

@inline function host_f_canon!(s)
    canonicalize!(s)
end

@inline function host_f_rref!(s)
    canonicalize_rref!(s)
end

@inline function device_f_canon!(
    s, synchronize;
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {block_SZ, batch_SZ}

    canonicalize!(
        s;
        block_size = block_size, batch_size = batch_size
        )
    synchronize()
end

@inline function device_f_rref!(
    s, synchronize;
    block_size::Val{block_SZ} = Val(default_block_size),
    batch_size::Val{batch_SZ} = Val(default_batch_size)
    ) where {block_SZ, batch_SZ}

    canonicalize_rref!(
        s;
        block_size = block_size, batch_size = batch_size
        )
    synchronize()
end

@inline function benchmark_KA_canonicalize(
    AT, synchronize, path;
    canon_type::Symbol = :canonicalize
    )

    host_time = zeros(Float64, length(n_qubits))
    device_time = zeros(Float64, length(batch_sizes), length(n_qubits))

    cache = AllocCache()
    for (i, n) in enumerate(n_qubits)
        @cached cache begin
            h_s = random_stabilizer(n)
            d_s = Stabilizer(
                AT(u32(h_s.phases)),
                h_s.nqubits,
                AT(reinterpret(UInt32, h_s.xzs))
                )
            synchronize()

            if canon_type == :canonicalize
                host_f = host_f_canon!
                device_f = device_f_canon!
            else
                host_f = host_f_rref!
                device_f = device_f_rref!
            end

            # Trigger compilation before benchmarking.
            host_f(copy(h_s))
            for (j, size) in enumerate(batch_sizes)
                device_f(copy(d_s), synchronize; batch_size = Val(size))
            end

            host_time[i] = @belapsed $host_f(copy($h_s)) evals = evals samples = samples seconds = seconds
            for (j, size) in enumerate(batch_sizes)
                device_time[j, i] =
                    @belapsed $device_f(copy($d_s), $synchronize; batch_size = Val($size)) evals = evals samples = samples seconds = seconds
            end
        end
    end
    unsafe_free!(cache)

    device_cat = [device_time[i, :] for i = 1 : length(batch_sizes)]
    title =
    "Performance uplift - $(canon_type)\n\
    Host threads = " * string(Threads.nthreads()) * " / " *
    string(Sys.CPU_THREADS) * ", Device block size = $default_block_size"
    xlabel = "Number of qubits"
    label = hcat(("Device - batch size = " .* string.(batch_sizes))..., "Host")

    plot(
        n_qubits, 10^3 .* hcat(device_cat..., host_time);
        shape = :circle, xticks = n_qubits, xscale = :log2, yscale = :log10,
        title = title, label = label, xlabel = xlabel, ylabel = "Runtime (ms)",
        background_color = :transparent
        )
    savefig("$path/runtime_$(canon_type).$format")

    plot(
        n_qubits, map(x -> host_time ./ x, device_cat);
        shape = :circle, xticks = n_qubits, xscale = :log2, title = title,
        label = hcat(label[1 : end - 1]...), xlabel = xlabel,
        ylabel = "Ratio (host/device)", background_color = :transparent
        )
    savefig("$path/ratio_$(canon_type).$format")
end

@inline function benchmark_KA_canonicalize_all(AT, synchronize, path)
    path = "$path/canonicalize"
    mkpath(path)
    benchmark_KA_canonicalize(AT, synchronize, path; canon_type = :canonicalize)
    benchmark_KA_canonicalize(AT, synchronize, path; canon_type = :rref)
end
