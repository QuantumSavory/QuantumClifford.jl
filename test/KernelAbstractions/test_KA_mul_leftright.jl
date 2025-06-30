using QuantumClifford: mul_left!, mul_right!, Tableau

function test_KA_mul_leftright(AT, synchronize)
	# Small sizes for encoding issues, large sizes for race conditions.
	test_sizes = [
		31, 32, 33, 63, 64, 65, 127, 128, 129,
		64 * 1023, 64 * 1024, 64 * 1025, 64 * 2047, 64 * 2048, 64 * 2049
		]

	for n in test_sizes
		for _ in 1:5
			h_p1 = random_pauli(n)
			d_p1 = PauliOperator(AT(h_p1.phase), h_p1.nqubits, AT(h_p1.xz))
			h_p2 = random_pauli(n)
			d_p2 = PauliOperator(AT(h_p2.phase), h_p2.nqubits, AT(h_p2.xz))
			h_s = Stabilizer(Tableau(
				rand(eltype(h_p1.phase), n) .& 0x3,
				n,
				rand(eltype(h_p1.xz), (length(h_p1.xz), n))
				))
			d_s = Stabilizer(Tableau(
				AT(h_s.tab.phases), h_s.tab.nqubits, AT(h_s.tab.xzs)
			))
			i = rand(1:n)

			d_o = mul_left!(copy(d_p2), d_p1)
			h_o = mul_left!(copy(h_p2), h_p1)
			synchronize()
			@test h_o.phase == Array(d_o.phase)
			@test h_o.xz == Array(d_o.xz)

			d_o = mul_right!(copy(d_p1), d_p2)
			h_o = mul_right!(copy(h_p1), h_p2)
			synchronize()
			@test h_o.phase == Array(d_o.phase)
			@test h_o.xz == Array(d_o.xz)

			d_L = mul_left!(copy(d_p2), d_p1)
			d_R = mul_right!(copy(d_p2), d_p1)
			synchronize()
			@test all(Array(d_L.phase) .== (-1)^comm(h_p1, h_p2) .* Array(d_R.phase))
			@test Array(d_L.xz) == Array(d_R.xz)

			d_o = mul_left!(copy(d_p2), d_s, i)
			h_o = mul_left!(copy(h_p2), h_s, i)
			synchronize()
			@test h_o.phase == Array(d_o.phase)
			@test h_o.xz == Array(d_o.xz)

			d_o = mul_right!(copy(d_p2), d_s, i)
			h_o = mul_right!(copy(h_p2), h_s, i)
			synchronize()
			@test h_o.phase == Array(d_o.phase)
			@test h_o.xz == Array(d_o.xz)

			d_o = mul_left!(copy(d_s), d_p2)
			h_o = mul_left!(copy(h_s), h_p2)
			synchronize()
			@test h_o.phases[i] == Array((@view d_o.phases[i]))
			@test h_o.xzs[:, i] == Array((@view d_o.xzs[:, i]))

			d_o = mul_right!(copy(d_s), d_p2)
			h_o = mul_right!(copy(h_s), h_p2)
			synchronize()
			@test h_o.phases[i] == Array((@view d_o.phases[i]))
			@test h_o.xzs[:, i] == Array((@view d_o.xzs[:, i]))

			# Potential race condition.
			d_o = mul_left!(copy(d_s), i, i)
			h_o = mul_left!(copy(h_s), i, i)
			synchronize()
			@test h_o.phases[i] == Array((@view d_o.phases[i]))
			@test h_o.xzs[:, i] == Array((@view d_o.xzs[:, i]))
		end
	end
end
