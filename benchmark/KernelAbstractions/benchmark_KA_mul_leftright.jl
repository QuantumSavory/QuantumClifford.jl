import QuantumClifford as QC
using GPUArrays: AllocCache, @cached, unsafe_free!
using BenchmarkTools: @belapsed
using Plots: scatter, savefig

@inline host_f(x, y; phases::Val{B} = Val(true)) where {B} =
	QC.mul_left!(x, y; phases = phases)

@inline function device_f(
	x, y, synchronize;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(QC.default_block_size),
	batch_size::Val{batch_SZ} = Val(QC.default_batch_size)
	) where {B, block_SZ, batch_SZ}

	QC.mul_left!(
		x, y;
		phases = phases, block_size = block_size, batch_size = batch_size
		)
	synchronize()

end

@inline function benchmark_KA_mul_leftright(
	AT, synchronize;
	phases::Val{B} = Val(true), platform_name = string(AT)
	) where {B}

	host_time = zeros(Float64, length(n_MiB))
	device_time = zeros(Float64, length(batch_sizes), length(n_MiB))

	# Keep the memory usage sane.
	cache = AllocCache()
	for (i, n) in enumerate(n_MiB)
		@cached cache begin
			# Each qubit requires 2 bits.
			h_p = QC.PauliOperator(
				zeros(Cuchar), n * MiB >> 1,
				zeros(UInt, cld(n * MiB, bit_count(UInt)))
				)
			d_p = QC.PauliOperator(AT(u32(h_p.phase)), h_p.nqubits, AT(h_p.xz))
			synchronize()
			# Trigger compilation before benchmarking.
			host_f(h_p, h_p; phases = phases)
			host_time[i] = @belapsed host_f($h_p, $h_p; phases = $phases)
			for (j, size) in enumerate(batch_sizes)
				device_f(
					d_p, d_p, synchronize;
					phases = phases, batch_size = Val(size)
					)
				device_time[j, i] =
					@belapsed device_f(
						$d_p, $d_p, $synchronize;
						phases = $phases, batch_size = Val($size)
						)
			end
		end
	end
	unsafe_free!(cache)

	device_cat = [device_time[i, :] for i = 1 : length(batch_sizes)]
	title = "Performance uplift (phases = $B)"
	xlabel = "Pauli operator size (MiB)"
	label = hcat(("Device - batch = " .* string.(batch_sizes))..., "Host")

	scatter(
		n_MiB, 10^3 .* hcat(device_cat..., host_time);
		xticks = n_MiB, xscale = :log2, yscale = :log10,
		title = title, label = label, xlabel = xlabel, ylabel = "Runtime (ms)"
		)
	savefig("runtime_" * platform_name * "_phase_$B.png")

	scatter(
		n_MiB, map(x -> host_time ./ x, device_cat);
		xticks = n_MiB, xscale = :log2, title = title,
		label = hcat(label[1 : end - 1]...), xlabel = xlabel,
		ylabel = "Ratio (host/device)"
		)
	savefig("ratio_" * platform_name * "_phase_$B.png")

end
