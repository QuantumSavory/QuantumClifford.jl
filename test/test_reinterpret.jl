@testitem "Reinterpret" begin
    using QuantumClifford

    # PauliOperator roundtrip via explicit constructor calls
    xbits = BitVector([true, true, false])
    zbits = BitVector([false, true, true])
    p = PauliOperator(xbits, zbits)
    p2 = reinterpret(UInt8, p)
    p3 = reinterpret(UInt64, p2)
    @test xbit(p) == xbit(p3)
    @test zbit(p) == zbit(p3)

    # Tableau roundtrip via explicit tableau construction
    pa1 = PauliOperator(BitVector([true, false]), BitVector([false, false])) # X X
    pa2 = PauliOperator(BitVector([false, false]), BitVector([true, true]))   # Z Z
    t = QuantumClifford.Tableau([pa1, pa2])
    t2 = reinterpret(UInt8, t)
    t3 = reinterpret(UInt64, t2)
    @test QuantumClifford.stab_to_gf2(t) == QuantumClifford.stab_to_gf2(t3)
end

@testitem "reinterpret edge cases" begin
    # After removing defensive layout checks, reinterpret should succeed and round-trip
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

@testset "reinterpret combinations" begin
    unsigned_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
    ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]

    for n in ns
        for Ti in unsigned_types, Tf in unsigned_types
            len = QuantumClifford._nchunks(n, Ti)
            xz = rand(Ti, len)
            p = PauliOperator(0x0, n, xz)
            try
                p2 = reinterpret(Tf, p)
                p3 = reinterpret(Ti, p2)
                @test xbit(p) == xbit(p3)
                @test zbit(p) == zbit(p3)
            catch e
                @test isa(e, ArgumentError)
            end
        end
    end
end


@testset "reinterpret extra edge cases" begin
    # 1) Zero-qubit PauliOperator (empty halves) should round-trip or at least raise ArgumentError
    p_zero = PauliOperator(BitVector(), BitVector())
    try
        pz2 = reinterpret(UInt64, p_zero)
        pz3 = reinterpret(eltype(p_zero.xz), pz2)
        @test xbit(p_zero) == xbit(pz3)
        @test zbit(p_zero) == zbit(pz3)
    catch e
        @test isa(e, ArgumentError)
    end

    # 2) Minimal backing: 1 byte total -> reinterpret to a 2-byte element should fail
    p_min = QuantumClifford.PauliOperator(0x0, 1, UInt8[0x0])
    @test_throws ArgumentError reinterpret(UInt16, p_min)

    # 3) Small backing -> reinterpret to a much larger element type should fail
    p_small2 = QuantumClifford.PauliOperator(0x0, 2, UInt8[0x0, 0x0])
    @test_throws ArgumentError reinterpret(UInt128, p_small2)

    # 4) Tableau with an odd number of rows (3 rows) -> reinterpret to a wider element should fail
    phases = UInt8[0x0]
    t_odd = QuantumClifford.Tableau(phases, 1, zeros(UInt8, 3, 1))
    @test_throws ArgumentError reinterpret(UInt16, t_odd)

    # 5) Successful tableau roundtrip for a simple shape (may throw ArgumentError depending on sizes)
    nch = QuantumClifford._nchunks(5, UInt8)
    xt = rand(UInt8, nch, 2)
    phases = zeros(UInt8, 2)
    t_ok = QuantumClifford.Tableau(phases, 5, xt)
    try
        t2 = reinterpret(UInt16, t_ok)
        t3 = reinterpret(eltype(t_ok.xzs), t2)
        @test QuantumClifford.stab_to_gf2(t_ok) == QuantumClifford.stab_to_gf2(t3)
    catch e
        @test isa(e, ArgumentError)
    end
end

@testset "tableau combinations" begin
    unsigned_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
    ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]
    rows_choices = (1, 2, 3)

    for n in ns
        for Ti in unsigned_types, Tf in unsigned_types
            for r in rows_choices
                nch = QuantumClifford._nchunks(n, Ti)
                xzs = rand(Ti, nch, r)
                phases = zeros(UInt8, r)
                t = QuantumClifford.Tableau(phases, n, xzs)
                try
                    t2 = reinterpret(Tf, t)
                    t3 = reinterpret(Ti, t2)
                    @test QuantumClifford.stab_to_gf2(t) == QuantumClifford.stab_to_gf2(t3)
                catch e
                    @test isa(e, ArgumentError)
                    msg = sprint(showerror, e)
                    @test occursin("cannot reinterpret", msg) && (occursin("backing bytes not divisible", msg) || occursin("resulting backing array length", msg))
                end
            end
        end
    end
end


@testset "stabilizer combinations" begin
    unsigned_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
    ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]
    rows_choices = (1, 2, 3)

    for n in ns
        for Ti in unsigned_types, Tf in unsigned_types
            for r in rows_choices
                nch = QuantumClifford._nchunks(n, Ti)
                xzs = rand(Ti, nch, r)
                phases = zeros(UInt8, r)
                t = QuantumClifford.Tableau(phases, n, xzs)
                s = QuantumClifford.Stabilizer(t)
                try
                    s2 = reinterpret(Tf, s)
                    s3 = reinterpret(Ti, s2)
                    @test QuantumClifford.stab_to_gf2(tab(s)) == QuantumClifford.stab_to_gf2(tab(s3))
                catch e
                    @test isa(e, ArgumentError)
                    msg = sprint(showerror, e)
                    @test occursin("cannot reinterpret", msg) && (occursin("backing bytes not divisible", msg) || occursin("resulting backing array length", msg))
                end
            end
        end
    end
end


@testset "pauliframe combinations" begin
    unsigned_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
    ns = [7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129]
    rows_choices = (1, 2, 3)

    for n in ns
        for Ti in unsigned_types, Tf in unsigned_types
            for r in rows_choices
                nch = QuantumClifford._nchunks(n, Ti)
                xzs = rand(Ti, nch, r)
                phases = zeros(UInt8, r)
                t = QuantumClifford.Tableau(phases, n, xzs)
                s = QuantumClifford.Stabilizer(t)
                pf = QuantumClifford.PauliFrame(s, falses(length(s), 2))
                try
                    pf2 = reinterpret(Tf, pf)
                    pf3 = reinterpret(Ti, pf2)
                    @test QuantumClifford.stab_to_gf2(tab(pf.frame)) == QuantumClifford.stab_to_gf2(tab(pf3.frame))
                    @test pf.measurements == pf3.measurements
                catch e
                    @test isa(e, ArgumentError)
                    msg = sprint(showerror, e)
                    @test occursin("cannot reinterpret", msg) && (occursin("backing bytes not divisible", msg) || occursin("resulting backing array length", msg))
                end
            end
        end
    end
end
