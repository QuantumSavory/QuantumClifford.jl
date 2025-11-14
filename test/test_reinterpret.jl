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

using Random

@testset "reinterpret edge cases" begin
    # After removing defensive layout checks, reinterpret should succeed and round-trip
    p_edge = P"XXXXX"
    p_edge2 = reinterpret(UInt128, p_edge)
    p_edge3 = reinterpret(eltype(p_edge.xz), p_edge2)
    @test xbit(p_edge) == xbit(p_edge3)
    @test zbit(p_edge) == zbit(p_edge3)

    small_xz = UInt8[0x0, 0x0, 0x0, 0x0]
    p_small = QuantumClifford.PauliOperator(0x0, 2, small_xz)
    p_small2 = reinterpret(UInt32, p_small)
    p_small3 = reinterpret(eltype(p_small.xz), p_small2)
    @test xbit(p_small) == xbit(p_small3)
    @test zbit(p_small) == zbit(p_small3)
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
                msg = sprint(showerror, e)
                allowed = (
                    "backing bytes not divisible by sizeof",
                    "backing bytes not divisible",
                    "is not even",
                    "insufficient to represent",
                )
                @test any(occursin(sub, msg) for sub in allowed)
            end
        end
    end
end

@testset "reinterpret failing edge cases" begin
    # create three paulis so the tableau has 3 rows -> total_bytes = 3 * sizeof(eltype)
    pa1 = PauliOperator(BitVector([true, false]), BitVector([false, false]))
    pa2 = PauliOperator(BitVector([false, false]), BitVector([true, true]))
    pa3 = PauliOperator(BitVector([true, true]), BitVector([true, false]))

    # PauliOperator with odd number of bytes (3 bytes) -- reinterpret to UInt16 should fail
    p_bad = QuantumClifford.PauliOperator(0x0, 1, UInt8[0x0, 0x0, 0x0])
    try
        reinterpret(UInt16, p_bad)
        @test false
    catch e
        @test isa(e, ArgumentError)
        msg = sprint(showerror, e)
        @test occursin("backing bytes not divisible", msg)
    end

    # Tableau backed by UInt8 with 3 rows -> total_bytes = 3 -> reinterpret to UInt16 should fail
    phases = UInt8[0x0]
    nqubits_tbl = 1
    xzs_bad = zeros(UInt8, 3, 1) # 3 rows * 1 byte each = 3 bytes -> not divisible by 2
    t_bad = QuantumClifford.Tableau(phases, nqubits_tbl, xzs_bad)
    try
        reinterpret(UInt16, t_bad)
        @test false
    catch e
        @test isa(e, ArgumentError)
        msg = sprint(showerror, e)
        @test occursin("backing bytes not divisible", msg)
    end

    # Stabilizer delegates to tableau reinterpret; should also fail when built from the bad tableau
    s_bad = QuantumClifford.Stabilizer(t_bad)
    try
        reinterpret(UInt16, s_bad)
        @test false
    catch e
        @test isa(e, ArgumentError)
        msg = sprint(showerror, e)
        @test occursin("backing bytes not divisible", msg)
    end

    # PauliFrame wraps a stabilizer/tableau; create a small measurements matrix and expect same failure
    pf_bad = QuantumClifford.PauliFrame(s_bad, falses(length(s_bad), 1))
    try
        reinterpret(UInt16, pf_bad)
        @test false
    catch e
        @test isa(e, ArgumentError)
        msg = sprint(showerror, e)
        @test occursin("backing bytes not divisible", msg)
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
                    @test occursin("backing bytes not divisible", msg) || occursin("is not even", msg) || occursin("insufficient to represent", msg)
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
                    @test occursin("backing bytes not divisible", msg) || occursin("is not even", msg) || occursin("insufficient to represent", msg)
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
                    @test occursin("backing bytes not divisible", msg) || occursin("is not even", msg) || occursin("insufficient to represent", msg)
                end
            end
        end
    end
end
