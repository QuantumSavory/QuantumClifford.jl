import QuantumClifford as QC, Plots
using BenchmarkTools: @belapsed

host_f(x, y; phases::Val{B} = Val(true)) where {B} =
	QC.mul_left!(x, y; phases = phases)

function device_f(
	l_phase, l_xz, r_phase, r_xz, synchronize;
	phases::Val{B} = Val(true),
	block_size::Val{block_SZ} = Val(QC.default_block_size),
	batch_size::Val{batch_SZ} = Val(QC.default_batch_size)
	) where {B, block_SZ, batch_SZ}

	QC.mul_device!(
		l_phase, l_xz, r_phase, r_xz;
		order_right_left = Val(true), phases = phases,
		block_size = block_size, batch_size = batch_size
		)
	synchronize()

end

function benchmark_KA_mul_leftright(
	AT, synchronize;
	phases::Val{B} = Val(true)
	platform_name = string(AT)
	) where {B}

	# There are 8 qubits in each byte.
	MiB = 8 * 1024 * 1024
	# Avoid consuming too much memory, 1 GiB is plenty.
	n_MiB = [2^i for i = 1:10]
	batch_sizes = [1, 4, 8, 16, 32, 64]
	host_time = zeros(Float64, length(n_MiB))
	device_time = zeros(Float64, length(batch_sizes), length(n_MiB))

	for (i, n) in enumerate(n_MiB)
		h_p = QC.random_pauli(n * MiB)
		d_p_xz = AT(h_p.xz)
		d_p_phase = AT(zeros(UInt32))
		synchronize()
		host_time[i] = @belapsed host_f($h_p, $h_p; phases = $phases)
		for (j, size) in enumerate(batch_sizes)
			device_f(
				d_p_phase, d_p_xz, d_p_phase, d_p_xz, synchronize;
				phases = phases, batch_size = Val(size)
				)
			device_time[j, i] =
				@belapsed device_f(
					$d_p_phase, $d_p_xz, $d_p_phase, $d_p_xz, $synchronize;
					phases = $phases, batch_size = Val($size)
					)
		end
		# Free up memory for the next round.
		h_p = nothing
		d_p_xz = nothing
		GC.gc()
	end

	device_cat = [device_time[i, :] for i = 1 : length(batch_sizes)]
	title = "Performance uplift (phases = $B)"
	xlabel = "Pauli operator size (MiB)"
	label = hcat(("Device - batch = " .* string.(batch_sizes))..., "Host")

	Plots.scatter(
		n_MiB, 10^3 .* hcat(device_cat..., host_time), xticks = n_MiB,
		xscale = :log2, yscale = :log10,
		title = title, label = label, xlabel = xlabel, ylabel = "Runtime (ms)"
		)
	Plots.savefig("runtime_" * platform_name * "_phase_$B.png")

	Plots.scatter(
		n_MiB, map(x -> host_time ./ x, device_cat), xticks = n_MiB,
		xscale = :log2, title = title, label = hcat(label[1 : end - 1]...),
		xlabel = xlabel, ylabel = "Ratio (host/device)"
		)
	Plots.savefig("ratio_" * platform_name * "_phase_$B.png")

end
