using QuantumClifford

using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
using QuantumClifford: projectremoverand!

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

@testset "Projective measurements" begin
    @testset "Stabilizer representation" begin
        s = S"XXX
                ZZI
                IZZ"
        ps, anticom, res = project!(copy(s), P"ZII")
        ps = canonicalize!(ps)
        @test anticom==1 && isnothing(res) && ps == S"ZII
                                                        IZI
                                                        IIZ"
        @test stab_looks_good(ps)

        ps, anticom, res = project!(copy(s), P"-XXX")
        @test anticom==0 && res[]==0x2 && ps == canonicalize!(copy(s))
        @test stab_looks_good(ps)

        ps, anticom, res = project!(copy(s), P"-XXX"; keep_result=false)
        @test anticom==0 && isnothing(res) && ps == s
        @test stab_looks_good(ps)

        for n in test_sizes
            s = random_stabilizer(n)
            m = random_pauli(n;nophase=true)
            mop = PauliMeasurement(m)
            ps, anticom, res = project!(copy(s),m)
            @test anticom==0x0 || ps[anticom]==m
            @test tab(ps).xzs == tab(apply!(copy(s),mop)).xzs
            @test stab_looks_good(ps)
            m = single_z(n,1)
            ps, anticom, res = project!(copy(s),m)
            @test anticom==0x0 || ps[anticom]==m
            @test stab_looks_good(ps)
            m = single_x(n,1)
            ps, anticom, res = project!(copy(s),m)
            @test anticom==0x0 || ps[anticom]==m
            @test stab_looks_good(ps)
        end
    end
    @testset "Destabilizer representation" begin
        for n in test_sizes
            s = canonicalize!(random_stabilizer(n))
            m = random_pauli(n;nophase=true)
            ps, anticom, res = project!(copy(s),m)
            dps, danticom, dres = project!(Destabilizer(copy(s)),m)
            @test destab_looks_good(dps)
            @test res==dres && canonicalize!(ps)==canonicalize!(stabilizerview(dps))
            m = single_z(n,1)
            ps, anticom, res = project!(copy(s),m)
            dps, danticom, dres = project!(Destabilizer(copy(s)),m)
            @test destab_looks_good(dps)
            @test res==dres && canonicalize!(ps)==canonicalize!(stabilizerview(dps))
            m = single_x(n,1)
            ps, anticom, res = project!(copy(s),m)
            dps, danticom, dres = project!(Destabilizer(copy(s)),m)
            @test destab_looks_good(dps)
            @test res==dres && canonicalize!(ps)==canonicalize!(stabilizerview(dps))
        end
    end
    @testset "Anticommutation indices and NA results" begin
        s = S" XXX
                -ZZI"
        ds = Destabilizer(copy(s))
        ms = MixedStabilizer(copy(s))
        mds = MixedDestabilizer(copy(s))

        p = P"IZZ"
        ps, a, r = project!(copy(s),p)
        @test stab_looks_good(ps)
        @test a==3 && isnothing(r)
        @test_throws BadDataStructure pds, a, r = project!(copy(ds),p)
        pms, a, r = project!(copy(ms),p)
        @test mixed_stab_looks_good(pms)
        @test pms.rank==3
        @test a==pms.rank && isnothing(r)
        pmds, a, r = project!(copy(mds),p)
        @test mixed_destab_looks_good(pmds)
        @test pmds.rank==3
        @test a==pmds.rank && isnothing(r)

        p = P"ZZI"
        ps, a, r = project!(copy(s),p)
        @test stab_looks_good(ps)
        @test a==0 && r==0x2
        @test_throws BadDataStructure pds, a, r = project!(copy(ds),p)
        pms, a, r = project!(copy(ms),p)
        @test mixed_stab_looks_good(pms)
        @test pms.rank==2
        @test a==0 && r==0x2
        pmds, a, r = project!(copy(mds),p)
        @test mixed_destab_looks_good(pmds)
        @test pmds.rank==2
        @test a==0 && r==0x2
        @test canonicalize!(ps)==canonicalize!(stabilizerview(pms))==canonicalize!(stabilizerview(pmds))

        p = P"XZZ"
        ps, a, r = project!(copy(s),p)
        @test stab_looks_good(ps)
        @test a==2 && isnothing(r)
        pds, a, r = project!(copy(ds),p)
        @test destab_looks_good(pds)
        @test a==2 && isnothing(r)
        pms, a, r = project!(copy(ms),p)
        @test mixed_stab_looks_good(pms)
        @test pms.rank==2
        @test a==2 && isnothing(r)
        pmds, a, r = project!(copy(mds),p)
        @test mixed_destab_looks_good(pmds)
        @test pmds.rank==2
        @test a==2 && isnothing(r)
        @test canonicalize!(ps)==canonicalize!(stabilizerview(pms))==canonicalize!(stabilizerview(pds))==canonicalize!(stabilizerview(pmds))
    end
    @testset "Mixed Destabilizer projection on logical operator" begin
        stab = one(MixedDestabilizer, 2,4)
        projzl = single_z(4,1)
        projzr = single_z(4,4)
        projxl = single_x(4,1)
        projxr = single_x(4,4)
        s, a, r = project!(copy(stab), projzl)
        @test mixed_destab_looks_good(s)
        @test a==0 && r==0x0       && stabilizerview(s)==S"Z___
                                                            _Z__"
        s, a, r = project!(copy(stab), projxl)
        @test mixed_destab_looks_good(s)
        @test a==1 && isnothing(r) && stabilizerview(s)==S"X___
                                                            _Z__"
        s, a, r = project!(copy(stab), projzr)
        @test mixed_destab_looks_good(s)
        @test a==3 && isnothing(r) && stabilizerview(s)==S"Z___
                                                            _Z__
                                                            ___Z"
        s, a, r = project!(copy(stab), projxr)
        @test mixed_destab_looks_good(s)
        @test a==3 && isnothing(r) && stabilizerview(s)==S"Z___
                                                            _Z__
                                                            ___X"
    end
    @testset "Interface Particularities" begin
        s = S"ZII IZI"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"IZI"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        s = S"ZII IZI III"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"IZI"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        s = MixedStabilizer(s, 2)
        ms, a, r = project!(copy(s), P"IZI")
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        @test ms.rank == 2
        ms, a, r = project!(copy(s), P"IIZ")
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        @test ms.rank == 3
        s = S"ZII IZI"
        s = Destabilizer(s)
        @test_throws BadDataStructure project!(copy(s), P"IZI"; keep_result=true)  # on comm
        @test_throws BadDataStructure project!(copy(s), P"IIZ"; keep_result=true)  # operators
        @test_throws BadDataStructure project!(copy(s), P"IZI"; keep_result=false) # in or out of
        @test_throws BadDataStructure project!(copy(s), P"IIZ"; keep_result=false) # the stabilizer
        s = S"ZII IZI IIZ"
        s = Destabilizer(s)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, 0x0)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing)
        s = S"ZII IZI"
        s = MixedDestabilizer(s)
        mds, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        @test mds.rank == 2
        mds, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        @test mds.rank == 3
        mds, a, r = project!(copy(s), P"IZI"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator in the stabilizer
        @test mds.rank == 2
        mds, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        @test mds.rank == 3
    end
    @testset "Results from canonicalization vs from destabilizer" begin
        @test generate!(P"_Z", S"XZ") === nothing # for bug fixed in 4b536231c3ee4e6446262fcc61ba8da669415bc8
        for n in test_sizes
            for r in [n, rand(n÷3:n*2÷3)]
                s = random_stabilizer(r,n)
                ms = MixedStabilizer(copy(s))
                d = Destabilizer(copy(s))
                md = MixedDestabilizer(copy(s))
                p = random_pauli(n,realphase=true)
                _, as, rs = project!(s,p)
                _, ams, rms = project!(ms,p)
                _, amd, rmd = project!(md,p)
                @test rs == rms == rmd
                @test (md.rank!=r) || (canonicalize!(s) == canonicalize!(stabilizerview(ms)))
                @test canonicalize!(stabilizerview(ms)) == canonicalize!(stabilizerview(md))
                if as == 0
                    @test ams == amd == 0
                end
                if r == n
                    _, ad, rd = project!(d,p)
                    @test s == canonicalize!(stabilizerview(d))
                    @test rd == rs
                    if as == 0
                        @test ad == 0
                    end
                end
            end
        end
    end
    @testset "Reported phase" begin
        s = S"ZII IZI"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer

        s = S"ZII -IZI"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer

        s = S"ZII IZI III"
        s = MixedStabilizer(s, 2)
        _, a, r = project!(copy(s), P"IZI")
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI")
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ")
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"-IIZ")
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        s = S"ZII -IZI III"
        s = MixedStabilizer(s, 2)
        _, a, r = project!(copy(s), P"IZI")
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI")
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer

        s = S"ZII IZI IIZ"
        s = Destabilizer(s)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, 0x0)
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (0, 0x2)
        s = S"ZII IZI -IIZ"
        s = Destabilizer(s)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, 0x2)
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (0, 0x0)

        s = S"ZII IZI"
        s = MixedDestabilizer(s)
        mds, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        mds, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        mds, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        mds, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (3, nothing) # on commuting operator out of the stabilizer
        s = S" ZII -IZI"
        s = MixedDestabilizer(s)
        mds, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        mds, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
    end
    @testset "Single Qubit Projections" begin
        for n in test_sizes
            s = MixedDestabilizer(random_destabilizer(n),n÷2);
            r = rand(1:n)
            px = single_x(n,r);
            py = single_y(n,r);
            pz = single_z(n,r);
            @test project!(copy(s),px) == projectX!(copy(s),r)
            @test project!(copy(s),py) == projectY!(copy(s),r)
            @test project!(copy(s),pz) == projectZ!(copy(s),r)
            sx = project!(copy(s),px)[1]
            sy = project!(copy(s),py)[1]
            sz = project!(copy(s),pz)[1]
            ssx = project!(copy(s),sMX(r))[1]
            ssy = project!(copy(s),sMY(r))[1]
            ssz = project!(copy(s),sMZ(r))[1]
            rssx = projectrand!(copy(s),sMX(r))[1]
            rssy = projectrand!(copy(s),sMY(r))[1]
            rssz = projectrand!(copy(s),sMZ(r))[1]
            @test project!(copy(sx),px) == projectX!(copy(sx),r)
            @test project!(copy(sy),py) == projectY!(copy(sy),r)
            @test project!(copy(sz),pz) == projectZ!(copy(sz),r)
            @test tab(apply!(copy(s),sMX(r))).xzs == tab(sx).xzs == tab(ssx).xzs == tab(rssx).xzs
            @test tab(apply!(copy(s),sMY(r))).xzs == tab(sy).xzs == tab(ssy).xzs == tab(rssy).xzs
            @test tab(apply!(copy(s),sMZ(r))).xzs == tab(sz).xzs == tab(ssz).xzs == tab(rssz).xzs
        end
    end
    @testset "projectremoverand!" begin
        for n in test_sizes
            n > 3 || continue
            state = MixedDestabilizer(random_destabilizer(n),n÷2)
            state_a = projectremoverand!(copy(state),projectX!,n)[1]
            state_b = traceout!(projectX!(state,n)[1],n)
            @test stab_to_gf2(canonicalize!(stabilizerview(state_a))) == stab_to_gf2(canonicalize!(stabilizerview(state_b)[:,1:end-1]))
        end
    end
    @testset "Redundant row permutations in `project!(::MixedDestabilizer)`" begin
        # Fixed in 41ed1d3c
        destab =  T"+ ZX_Y_YXZ
                    + XY_Y____
                    + _Z_XXY__
                    + _ZYXXY__
                    + X__Y_ZXZ
                    + X__YXZXZ
                    + ___YXXZZ
                    + _______Z"
        stab =    T"+ X_______
                    + _X_Y____
                    + __ZY____
                    + __Z_____
                    + ___YZY__
                    + X__YZYZZ
                    + X____YZZ
                    + ______YX"
        t = MixedDestabilizer(vcat(destab,stab), 8)
        @test mixed_destab_looks_good(t)
        c = copy(stabilizerview(t)[[1,3,5,7]])
        traceout!(t,[1,4,3,6])
        @test mixed_destab_looks_good(t)
        project!(t,c[1])
        @test mixed_destab_looks_good(t)
        project!(t,c[2])
        @test mixed_destab_looks_good(t) # This used to fail because anticomlog==rank+1 leading to a repeated row permutation
    end
end
