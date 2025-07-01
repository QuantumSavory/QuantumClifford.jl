import QuantumClifford as QC
using GPUArrays: AllocCache, @cached, unsafe_free!

@inline function test_KA_mul_leftright(AT, synchronize)
	# Small sizes for encoding issues, large sizes for race conditions.
	test_sizes = [
		31, 32, 33, 63, 64, 65, 127, 128, 129,
		64 * 1023, 64 * 1024, 64 * 1025, 64 * 2047, 64 * 2048, 64 * 2049
		]
	u32(z) = map(x -> UInt32(x), z)

	cache = AllocCache()
	for n in test_sizes
		# Keep the memory usage sane.
		rows = min(1024, n)
		for _ in 1:16
			@cached cache begin
				h_p1 = QC.random_pauli(n)
				d_p1 = QC.PauliOperator(
					AT(u32(h_p1.phase)), h_p1.nqubits, AT(h_p1.xz)
					)
				h_p2 = QC.random_pauli(n)
				d_p2 = QC.PauliOperator(
					AT(u32(h_p2.phase)), h_p2.nqubits, AT(h_p2.xz)
					)
				h_s = QC.Stabilizer(
					QC.Tableau(
						rand(eltype(h_p1.phase), rows) .& 0x3,
						n,
						rand(eltype(h_p1.xz), (length(h_p1.xz), rows))
						)
					)
				d_s = QC.Stabilizer(
					QC.Tableau(
						AT(u32(h_s.tab.phases)),
						h_s.tab.nqubits,
						AT(h_s.tab.xzs)
						)
					)
				i = rand(1:rows)

				d_o = QC.mul_left!(copy(d_p2), d_p1)
				h_o = QC.mul_left!(copy(h_p2), h_p1)
				synchronize()
				@test begin
					h_o.phase == Array(d_o.phase)
					h_o.xz == Array(d_o.xz)
				end

				d_o = QC.mul_right!(copy(d_p1), d_p2)
				h_o = QC.mul_right!(copy(h_p1), h_p2)
				synchronize()
				@test begin
					h_o.phase == Array(d_o.phase)
					h_o.xz == Array(d_o.xz)
				end

				d_L = QC.mul_left!(copy(d_p2), d_p1)
				d_R = QC.mul_right!(copy(d_p2), d_p1)
				synchronize()
				@test begin
					all(
						(Array(d_L.phase) .- Array(d_R.phase)) .& 0x3
						.== 2 * QC.comm(h_p1, h_p2)
						)
					Array(d_L.xz) == Array(d_R.xz)
				end

				d_o = QC.mul_left!(copy(d_p2), d_s, i)
				h_o = QC.mul_left!(copy(h_p2), h_s, i)
				synchronize()
				@test begin
					h_o.phase == Array(d_o.phase)
					h_o.xz == Array(d_o.xz)
				end

				d_o = QC.mul_right!(copy(d_p2), d_s, i)
				h_o = QC.mul_right!(copy(h_p2), h_s, i)
				synchronize()
				@test begin
					h_o.phase == Array(d_o.phase)
					h_o.xz == Array(d_o.xz)
				end

				d_o = QC.mul_left!(copy(d_s), d_p2)
				h_o = QC.mul_left!(copy(h_s), h_p2)
				synchronize()
				@test begin
					(@view h_o.tab.phases[i]) ==
						Array((@view d_o.tab.phases[i]))
					(@view h_o.tab.xzs[:, i]) ==
						Array((@view d_o.tab.xzs[:, i]))
				end

				d_o = QC.mul_right!(copy(d_s), d_p2)
				h_o = QC.mul_right!(copy(h_s), h_p2)
				synchronize()
				@test begin
					(@view h_o.tab.phases[i]) ==
						Array((@view d_o.tab.phases[i]))
					(@view h_o.tab.xzs[:, i]) ==
						Array((@view d_o.tab.xzs[:, i]))
				end

				# Potential race condition.
				d_o = QC.mul_left!(copy(d_s), i, i)
				h_o = QC.mul_left!(copy(h_s), i, i)
				synchronize()
				@test begin
					(@view h_o.tab.phases[i]) ==
						Array((@view d_o.tab.phases[i]))
					(@view h_o.tab.xzs[:, i]) ==
						Array((@view d_o.tab.xzs[:, i]))
				end
			end
		end
	end
	unsafe_free!(cache)
end
