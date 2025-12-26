@testitem "Reinterpret" begin
	using QuantumClifford

	# PauliOperator with explicit UInt64 storage (8 bytes per element)
	# For 3 qubits: xbits = [true, true, false], zbits = [false, true, true]
	# Encoded as: X bits in first chunk (0b011 = 3), Z bits in second chunk (0b110 = 6)
	p = PauliOperator(0x0, 3, UInt64[3, 6])
	@test xbit(p) == [true, true, false]
	@test zbit(p) == [false, true, true]

	# Tableau with explicit UInt64 storage
	# For 2 qubits: row1 = X₁X₂ (xbits=0b11=3, zbits=0b00=0), row2 = Z₁Z₂ (xbits=0b00=0, zbits=0b11=3)
	nch = QuantumClifford._nchunks(2, UInt64)
	xzs = reshape(UInt64[3, 0, 0, 3], nch, 2)
	phases = UInt8[0x0, 0x0]
	t = QuantumClifford.Tableau(phases, 2, xzs)
	@test QuantumClifford.stab_to_gf2(t)[1, :] == [true, true, false, false]  # X₁X₂
	@test QuantumClifford.stab_to_gf2(t)[2, :] == [false, false, true, true]  # Z₁Z₂
end

@testitem "reinterpret edge cases" begin
	p_edge = P"XXXXX"
	try
		p_edge2 = reinterpret(UInt128, p_edge)
		p_edge3 = reinterpret(eltype(p_edge.xz), p_edge2)
		@test xbit(p_edge) == xbit(p_edge3)
		@test zbit(p_edge) == zbit(p_edge3)
	catch e
		@test isa(e, ArgumentError)
	end

	small_xz = UInt8[0x0, 0x0, 0x0, 0x0]
	p_small = QuantumClifford.PauliOperator(0x0, 2, small_xz)
	try
		p_small2 = reinterpret(UInt32, p_small)
		p_small3 = reinterpret(eltype(p_small.xz), p_small2)
		@test xbit(p_small) == xbit(p_small3)
		@test zbit(p_small) == zbit(p_small3)
	catch e
		@test isa(e, ArgumentError)
	end
end
	@test t[1, :] == [true, true, false, false]  # X₁X₂
	@test t[2, :] == [false, false, true, true]  # Z₁Z₂
	unsigned_types = subtypes(Unsigned)
	ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]

	for n in ns
		for Ti in unsigned_types, Tf in unsigned_types
			p = zero(PauliOperator, n)
			random_pauli!(p)
			try
				p2 = reinterpret(Tf, p)
				p3 = reinterpret(Ti, p2)
				@test xbit(p) == xbit(p3)
				@test zbit(p) == zbit(p3)
			catch e
				@test xbit(p) == xbit(p3)
				@test zbit(p) == zbit(p3)
		end
	end
end


@testset "reinterpret extra edge cases" begin
	p_zero = PauliOperator(BitVector(), BitVector())
	try
		pz2 = reinterpret(UInt64, p_zero)
		pz3 = reinterpret(eltype(p_zero.xz), pz2)
		@test xbit(p_zero) == xbit(pz3)
		@test zbit(p_zero) == zbit(pz3)
	catch e
		@test isa(e, ArgumentError)
	end

	p_min = QuantumClifford.PauliOperator(0x0, 1, UInt8[0x0])
				@test t == t3

	p_small2 = QuantumClifford.PauliOperator(0x0, 2, UInt8[0x0, 0x0])
	@test_throws ArgumentError reinterpret(UInt128, p_small2)

	phases = UInt8[0x0]
	t_odd = QuantumClifford.Tableau(phases, 1, zeros(UInt8, 3, 1))
	@test_throws ArgumentError reinterpret(UInt16, t_odd)

	nch = QuantumClifford._nchunks(5, UInt8)
	xt = rand(UInt8, nch, 2)
	phases = zeros(UInt8, 2)
	t_ok = QuantumClifford.Tableau(phases, 5, xt)
	try
		t2 = reinterpret(UInt16, t_ok)
		t3 = reinterpret(eltype(t_ok.xzs), t2)
		@test t_ok == t3
				@test tab(s) == tab(s3)
		@test isa(e, ArgumentError)
	end
end

@testset "tableau combinations" begin
	unsigned_types = subtypes(Unsigned)
	ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]
	rows_choices = (1, 2, 3)

	for n in ns
		for Ti in unsigned_types, Tf in unsigned_types
			for r in rows_choices
				t = zero(Tableau, r, n)
				random_tableau!(t)
				try
			@test tab(d) == tab(d3)
					t3 = reinterpret(Ti, t2)
					@test t == t3
				catch e
					@test isa(e, ArgumentError)
				end
			end
		end
	end
end


@testset "stabilizer combinations" begin
	unsigned_types = subtypes(Unsigned)
	ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]
	rows_choices = (1, 2, 3)
			@test tab(ms) == tab(ms3)
	for n in ns
		for Ti in unsigned_types, Tf in unsigned_types
			for r in rows_choices
				t = zero(Tableau, r, n)
				random_tableau!(t)
				s = QuantumClifford.Stabilizer(t)
				try
					s2 = reinterpret(Tf, s)
					s3 = reinterpret(Ti, s2)
					@test s == s3
				catch e
					@test isa(e, ArgumentError)
				end
			end
		end
			@test tab(md) == tab(md3)
end


@testset "destabilizer combinations" begin
	unsigned_types = subtypes(Unsigned)
	ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]

	for n in ns
		for Ti in unsigned_types, Tf in unsigned_types
			t = zero(Tableau, n, n)
			random_tableau!(t)
			s = QuantumClifford.Stabilizer(t)
			d = QuantumClifford.Destabilizer(s)
			try
				d2 = reinterpret(Tf, d)
				d3 = reinterpret(Ti, d2)
				@test tab(pf.frame) == tab(pf3.frame)
			catch e
				@test isa(e, ArgumentError)
			end
		end
	end
end


@testset "mixedstabilizer combinations" begin
	unsigned_types = subtypes(Unsigned)
	ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]

	for n in ns
		for Ti in unsigned_types, Tf in unsigned_types
			r = min(n-1, max(1, n÷2)) # Make it rank-deficient
			t = zero(Tableau, r, n)
			random_tableau!(t)
			s = QuantumClifford.Stabilizer(t)
			ms = QuantumClifford.MixedStabilizer(s)
			try
				ms2 = reinterpret(Tf, ms)
				ms3 = reinterpret(Ti, ms2)
				@test ms == ms3
				@test rank(ms) == rank(ms3)
			catch e
				@test isa(e, ArgumentError)
			end
		end
	end
end


@testset "mixeddestabilizer combinations" begin
	unsigned_types = subtypes(Unsigned)
	ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]

	for n in ns
		for Ti in unsigned_types, Tf in unsigned_types
			r = min(n-1, max(1, n÷2)) # Make it rank-deficient
			t = zero(Tableau, r, n)
			random_tableau!(t)
			s = QuantumClifford.Stabilizer(t)
			md = QuantumClifford.MixedDestabilizer(s)
			try
				md2 = reinterpret(Tf, md)
				md3 = reinterpret(Ti, md2)
				@test md == md3
				@test rank(md) == rank(md3)
			catch e
				@test isa(e, ArgumentError)
			end
		end
	end
end


@testset "pauliframe combinations" begin
	unsigned_types = subtypes(Unsigned)
	ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]
	rows_choices = (1, 2, 3)

	for n in ns
		for Ti in unsigned_types, Tf in unsigned_types
			for r in rows_choices
				t = zero(Tableau, r, n)
				random_tableau!(t)
				s = QuantumClifford.Stabilizer(t)
				pf = QuantumClifford.PauliFrame(s, falses(length(s), 2))
				try
					pf2 = reinterpret(Tf, pf)
					pf3 = reinterpret(Ti, pf2)
					@test pf == pf3
				catch e
					@test isa(e, ArgumentError)
				end
			end
		end
	end
end
