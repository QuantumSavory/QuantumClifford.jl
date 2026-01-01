@testitem "reinterpret" begin
    using QuantumClifford
    using Test
    using InteractiveUtils

    # Local aliases for brevity
    const Tableau = QuantumClifford.Tableau
    const random_tableau! = QuantumClifford.random_tableau!
    const random_pauli! = QuantumClifford.random_pauli!

    reinterpret_error_matches(e, needle="Unable to reinterpret") = begin
        @test (isa(e, ArgumentError) || isa(e, MethodError))
        needles = needle isa AbstractString ? (needle,) : needle
        @test any(n -> occursin(n, sprint(showerror, e)), needles)
    end

    @testset "basic pauli/tableau roundtrips" begin
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

    @testset "reinterpret edge cases" begin
        p_edge = P"XXXXX"
        try
            p_edge2 = reinterpret(UInt128, p_edge)
            p_edge3 = reinterpret(eltype(p_edge.xz), p_edge2)
            @test xbit(p_edge) == xbit(p_edge3)
            @test zbit(p_edge) == zbit(p_edge3)
        catch e
            reinterpret_error_matches(e, "Unable to reinterpret pauli storage")
        end

        small_xz = UInt8[0x0, 0x0, 0x0, 0x0]
        p_small = QuantumClifford.PauliOperator(0x0, 2, small_xz)
        try
            p_small2 = reinterpret(UInt32, p_small)
            p_small3 = reinterpret(eltype(p_small.xz), p_small2)
            @test xbit(p_small) == xbit(p_small3)
            @test zbit(p_small) == zbit(p_small3)
        catch e
            reinterpret_error_matches(e, "Unable to reinterpret pauli storage")
        end
    end

    @testset "pauli combinations" begin
        unsigned_types = [UInt8, UInt64, UInt128]
        ns = [64]

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
                    reinterpret_error_matches(e, "Unable to reinterpret pauli storage")
                end
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
            reinterpret_error_matches(e, "Unable to reinterpret pauli storage")
        end

        p_min = QuantumClifford.PauliOperator(0x0, 1, UInt8[0x0])
        @test_throws "Unable to reinterpret pauli storage" reinterpret(UInt8, p_min)

        p_small2 = QuantumClifford.PauliOperator(0x0, 2, UInt8[0x0, 0x0])
        @test_throws "Unable to reinterpret pauli storage" reinterpret(UInt128, p_small2)

        phases = UInt8[0x0]
        t_odd = QuantumClifford.Tableau(phases, 1, zeros(UInt8, 3, 1))
        @test_throws "Unable to reinterpret tableau storage" reinterpret(UInt16, t_odd)

        nch = QuantumClifford._nchunks(5, UInt8)
        xt = rand(UInt8, nch, 2)
        phases = zeros(UInt8, 2)
        t_ok = QuantumClifford.Tableau(phases, 5, xt)
        try
            t2 = reinterpret(UInt16, t_ok)
            t3 = reinterpret(eltype(t_ok.xzs), t2)
            @test t_ok == t3
        catch e
            reinterpret_error_matches(e, "Unable to reinterpret tableau storage")
        end
    end

    @testset "tableau combinations" begin
        unsigned_types = [UInt8, UInt64, UInt128]
        ns = [64]
        rows_choices = (3,)

        for n in ns
            for Ti in unsigned_types, Tf in unsigned_types
                for r in rows_choices
                    t = zero(Tableau, r, n)
                    random_tableau!(t)
                    try
                        t2 = reinterpret(Tf, t)
                        t3 = reinterpret(Ti, t2)
                        @test t == t3
                    catch e
                        reinterpret_error_matches(e, "Unable to reinterpret tableau storage")
                    end
                end
            end
        end
    end

    @testset "tableau layout variations" begin
        # Test that layout functions work with reinterpret
        Ti, Tf = UInt8, UInt64
        n, r = 7, 2
        for layout_fn in (identity, fastrow, fastcolumn)
            t = layout_fn(zero(Tableau, r, n))
            random_tableau!(t)
            try
                t2 = reinterpret(Tf, t)
                t3 = reinterpret(Ti, t2)
                @test t == t3
            catch e
                reinterpret_error_matches(e, "Unable to reinterpret tableau storage")
            end
        end
    end

    @testset "stabilizer combinations" begin
        unsigned_types = [UInt8, UInt64, UInt128]
        ns = [64]
        rows_choices = (3,)

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
                        reinterpret_error_matches(e, "Unable to reinterpret stabilizer storage")
                    end
                end
            end
        end
    end

    @testset "destabilizer combinations" begin
        unsigned_types = [UInt8, UInt64, UInt128]
        ns = [64]

        for n in ns
            for Ti in unsigned_types, Tf in unsigned_types
                t = zero(Tableau, n, n)
                random_tableau!(t)
                s = QuantumClifford.Stabilizer(t)
                d = QuantumClifford.Destabilizer(s)
                try
                    d2 = reinterpret(Tf, d)
                    d3 = reinterpret(Ti, d2)
                    @test tab(d) == tab(d3)
                catch e
                    reinterpret_error_matches(e, "Unable to reinterpret destabilizer storage")
                end
            end
        end
    end

    @testset "mixedstabilizer combinations" begin
        unsigned_types = [UInt8, UInt64, UInt128]
        ns = [64]

        for n in ns
            for Ti in unsigned_types, Tf in unsigned_types
                r = min(n-1, max(1, n÷2))
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
                    reinterpret_error_matches(e, "Unable to reinterpret mixedstabilizer storage")
                end
            end
        end
    end

    @testset "mixeddestabilizer combinations" begin
        unsigned_types = [UInt8, UInt64, UInt128]
        ns = [64]

        for n in ns
            for Ti in unsigned_types, Tf in unsigned_types
                r = min(n-1, max(1, n÷2))
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
                    reinterpret_error_matches(e, "Unable to reinterpret mixeddestabilizer storage")
                end
            end
        end
    end

    @testset "pauliframe combinations" begin
        unsigned_types = [UInt8, UInt64, UInt128]
        ns = [64]
        rows_choices = (3,)

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
                        reinterpret_error_matches(e, ("Unable to reinterpret pauliframe storage", "Unable to reinterpret tableau storage", "Cannot `convert`"))
                    end
                end
            end
        end
    end
end
