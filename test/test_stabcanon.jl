using QuantumClifford

using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

@testset "Stabilizer canonicalization" begin
    @testset "Default canonicalization" begin
        s = S"- XZZZZ_____
            - _YZY___YX_
            - __XXZ__YX_
            + Z_Z_Y__YXZ
            + _____Z____
            + __________
            + __________
            + ______YYX_
            + __________
            + __________"
        canonicalize!(s)
        t = S"- XZZZZ_____
            - _YZY___YX_
            - __XXZ__YX_
            + Z_Z_Y__YXZ
            + ______YYX_
            + _____Z____
            + __________
            + __________
            + __________
            + __________"
        @test s == t
        for n in test_sizes
            @test stab_looks_good(random_stabilizer(n))
        end
    end
    @testset "Gottesman canonicalization" begin
        for n in test_sizes
            for nrows in [n, rand(n÷3:n*2÷3)]
                rs = random_stabilizer(nrows,n)
                c = canonicalize!(copy(rs))
                g, _, _, perm1, perm2 = canonicalize_gott!(copy(rs))
                c1 = canonicalize!(permute!(permute!(copy(rs),perm1),perm2))
                cg = canonicalize!(copy(g))
                @test cg == c1
                @test stab_looks_good(g)
            end
        end
    end
    @testset "Canonicalization of complex tableaus" begin
        for n in test_sizes
            for nrows in [n, rand(n÷3:n*2÷3)]
                rs = random_stabilizer(nrows,n)
                rs_m = MixedStabilizer(copy(rs))
                rs_d = Destabilizer(copy(rs))
                rs_md = MixedDestabilizer(copy(rs))
                c  = canonicalize!(copy(rs))
                mc  = canonicalize!(copy(rs_m))
                dc = canonicalize!(copy(rs_d))
                mdc = canonicalize!(copy(rs_md))
                @test stabilizerview(mdc) == stabilizerview(mc) == stabilizerview(dc) == c
                @test stab_looks_good(c)
                @test mixed_stab_looks_good(mc)
                @test destab_looks_good(dc)
                @test mixed_destab_looks_good(mdc)
                c  = canonicalize_rref!(copy(rs),1:n)[1]
                mc  = canonicalize_rref!(copy(rs_m),1:n)[1]
                dc = canonicalize_rref!(copy(rs_d),1:n)[1]
                mdc = canonicalize_rref!(copy(rs_md),1:n)[1]
                @test stabilizerview(mdc) == stabilizerview(mc) == stabilizerview(dc) == c
                @test stab_looks_good(c)
                @test mixed_stab_looks_good(mc)
                @test destab_looks_good(dc)
                @test mixed_destab_looks_good(mdc)
            end
        end
    end
end
