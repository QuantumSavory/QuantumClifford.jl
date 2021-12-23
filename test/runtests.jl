using Test, Random, Documenter
using QuantumClifford
using QuantumClifford: stab_looks_good, mixed_stab_looks_good, destab_looks_good, mixed_destab_looks_good
using QuantumClifford: CNOTcol, SWAPcol, Hadamardcol, Phasecol, CliffordIdcol, mul_left!
using QuantumClifford.Experimental.NoisyCircuits
using Quantikz: circuit2string, QuantikzOp
using Nemo
using LinearAlgebra: inv
import AbstractAlgebra

test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

function doset(descr)
    if length(ARGS) == 0
        return true
    end
    for a in ARGS
        if occursin(lowercase(a), lowercase(descr))
            return true
        end
    end
    return false
end

function doctests()
if doset("Doctests") && VERSION >= v"1.6"
@testset "Doctests" begin
    DocMeta.setdocmeta!(QuantumClifford, :DocTestSetup, :(using QuantumClifford); recursive=true)
    doctest(QuantumClifford)
end
end
end

function test_paulistab()
if doset("Pauli Operators")
@testset "Pauli Operators" begin
    @testset "Parsing, constructors, and properties" begin
        @test P"-iXYZ" == PauliOperator(0x3, 3, vcat(BitArray([1,1,0]).chunks, BitArray([0,1,1]).chunks))
        @test P"-iXYZ" == PauliOperator(0x3, Bool[1,1,0], Bool[0,1,1])
        @test xbit(P"-iXYZ") == Bool[1,1,0]
        @test zbit(P"-iXYZ") == Bool[0,1,1]
        @test P"-iXYZ".xz == UInt64[0x03, 0x06]
        @test P"-iXYZ".phase[] == 0x03
        @test P"-iXYZ".nqubits == 3
        @test size(P"-iXYZ") == (3,)
    end
    @testset "Indexing" begin
        @test eachindex(P"IXYZ") == 1:4
        @test P"IXYZ"[3] == (true, true)
        p = P"IXYZ"
        @test p[[3,2]] == P"YX"
        p[4] = (true,false)
        @test p == P"IXYX"
    end
    @testset "Elementary operations" begin
        @test P"X"*P"Z" == P"-iY"
        @test comm(P"XX",P"YY") == 0x0
        @test comm(P"XZ",P"YZ") == 0x1
        @test prodphase(P"XX",P"YY") == 0x2
        @test prodphase(P"ZZZ",P"XXX") == 0x3
    end
    @testset "Commutation implies real phase" begin
        for i in 1:10
            for n in test_sizes
                p1,p2 = random_pauli(n; nophase=true), random_pauli(n; nophase=true)
                com = comm(p1,p2)==0x0
                p = prodphase(p1,p2)
                rea = p==0x0 || p==0x2
                @test (com && rea) || (!com && !rea)
            end
        end
    end
    @testset "Prodphase" begin
        for i in 1:10
            for n in test_sizes
                p1,p2 = random_pauli(n; nophase=true), random_pauli(n; nophase=true)
                p = prodphase(p1.xz,p2.xz)
                @test p == QuantumClifford._stim_prodphase(p1.xz,p2.xz)&0x3
            end
        end
    end
    @testset "Single qubit Paulis and their action" begin
        for i in 1:3
            for n in test_sizes
                for t in [Stabilizer, Destabilizer, MixedDestabilizer]
                    ix, iy, iz = rand(1:n), rand(1:n), rand(1:n)
                    px = single_x(n,ix)
                    py = single_y(n,iy)
                    pz = single_z(n,iz)
                    s1 = t(random_stabilizer(n))
                    s2 = copy(s1)
                    apply!(s1,px)
                    apply_single_x!(s2,ix)
                    @test s1==s2
                    apply!(s1,py)
                    apply_single_y!(s2,iy)
                    @test s1==s2
                    apply!(s1,pz)
                    apply_single_z!(s2,iz)
                    @test s1==s2
                end
            end
        end
    end
end
end

if doset("Pure and Mixed state initialization")
@testset "Pure and Mixed state initialization" begin
    @testset "Destabilizer initialization" begin
        for n in test_sizes
            @test destab_looks_good(Destabilizer(random_stabilizer(n)))
        end
    end
    @testset "Mixed destabilizer initialization" begin
        for n in test_sizes
            @test n<10 || mixed_destab_looks_good(MixedDestabilizer(random_stabilizer(rand(n÷2+1:n-4),n)))
        end
        # Test initialization out of overdetermined stabs
        stabs = [S"XX
                   II",
                 S"XX
                   XX",
                 S"ZZ
                   ZZ"]
        for s in stabs
            md = MixedDestabilizer(s[1:1])
            @test_broken MixedDestabilizer(s) == md
        end
    end
    @testset "Tensor products over stabilizers" begin
        for n in test_sizes
            n<10 && continue
            l = random_stabilizer(rand(n÷2+1:n-2),n)
            r = random_stabilizer(rand(n÷3:n÷2),rand(n÷2:n))
            s = l⊗r
            ds = MixedDestabilizer(l)⊗MixedDestabilizer(r)
            @test mixed_destab_looks_good(ds)
            canonicalize!(s)
            dss = canonicalize!(copy(stabilizerview(ds)))
            @test s == dss
        end
    end
end
end

if doset("Stabilizer canonicalization")
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
                c1 = canonicalize!(colpermute!(colpermute!(copy(rs),perm1),perm2))
                cg = canonicalize!(copy(g))
                @test cg == c1
                @test stab_looks_good(g)
            end
        end
    end
    @testset "Canonicalization of complex tableaus" begin
        for n in test_sizes
            for nrows in [n, rand(n÷3:n*2÷3)]
                if nrows==0
                    @test_broken error("can not process empty stab")
                    continue
                end
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
end

if doset("Stabilizer indexing")
@testset "Stabilizer indexing" begin
    s = random_stabilizer(10,10)
    @test s[1,1] == s[[1,3,4],[1,3,5]][1,1]
    @test s[1,1] == s[:,[1,3,5]][1,1]
    @test s[1,1] == s[1,[1,3,5]][1]
    @test s[1,1] == s[[1,3,5],:][1,1]
    @test s[1,1] == s[[1,3,5],1][1,1]
    @test s[1,1] == s[:,1][1,1]
    @test s[1,1] == s[1,:][1]
    @test s[1,1] == s[:,:][1,1]
    @test isa(s[1], PauliOperator)
    @test isa(s[1,:], PauliOperator)
    @test isa(s[1,[1,2,3]], PauliOperator)
end
end
if doset("Low-level tableaux ops")
@testset "Low-level tableaux ops" begin
    for n in test_sizes
        p1 = random_pauli(n)
        p11 = copy(p1)
        p2 = random_pauli(n)
        p3 = p2*p1
        s = Stabilizer([p1,p2])
        @test QuantumClifford._mul_left_nonvec!(copy(p1).xz,p2.xz)&3 == mul_left!(copy(p1).xz,p2.xz)&3 
        @test prodphase(p2,p1) == mul_left!(p1,p2).phase[]
        mul_left!(p11,s,2)
        mul_left!(s,1,2)
        @test p1 == p11 == p3 == s[1] 
    end
end
end
if doset("GF(2) representations")
@testset "GF(2) representations" begin
    @testset "Equivalence of GF(2) Gaussian elimination and Stabilizer canonicalization" begin
        for n in test_sizes
            for rep in 1:5
                s = random_stabilizer(n)[randperm(n)[1:rand(n÷2+1:n)]]
                cs = canonicalize!(copy(s));
                H = stab_to_gf2(cs);
                cH = gf2_gausselim!(stab_to_gf2(s));
                @test H==cH
            end
        end
    end
    @testset "GF(2) H and G matrices" begin
        for n in test_sizes
            for rep in 1:5
                H = random_invertible_gf2(n)[randperm(n)[1:rand(n÷2+1:n)],:]
                H = gf2_gausselim!(H)
                G = gf2_H_to_G(H)
                @test sum(G*H' .%2)==0;
            end
        end
    end
end
end
end

function test_operations()
if doset("Projective measurements")
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
            ps, anticom, res = project!(copy(s),m)
            @test anticom==0x0 || ps[anticom]==m
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
        @test a==0 && isnothing(r)
        @test_throws BadDataStructure pds, a, r = project!(copy(ds),p)
        pms, a, r = project!(copy(ms),p)
        @test mixed_stab_looks_good(pms)
        @test pms.rank==3
        @test a==0 && isnothing(r)
        pmds, a, r = project!(copy(mds),p)
        @test mixed_destab_looks_good(pmds)
        @test pmds.rank==3
        @test a==0 && isnothing(r)

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
        @test a==0 && isnothing(r) && stabilizerview(s)==S"Z___
                                                           _Z__
                                                           ___Z"
        s, a, r = project!(copy(stab), projxr)
        @test mixed_destab_looks_good(s)
        @test a==0 && isnothing(r) && stabilizerview(s)==S"Z___
                                                           _Z__
                                                           ___X"
    end
    @testset "Interface Particularities" begin
        s = S"ZII
              IZI"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"IZI"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        s = S"ZII
              IZI
              III"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"IZI"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        s = MixedStabilizer(s, 2)
        ms, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        @test ms.rank == 2
        ms, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        @test ms.rank == 3
        ms, a, r = project!(copy(s), P"IZI"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator in the stabilizer
        @test ms.rank == 2
        ms, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        @test ms.rank == 3
        s = S"ZII
              IZI"
        s = Destabilizer(s)
        @test_throws BadDataStructure project!(copy(s), P"IZI"; keep_result=true)  # on comm
        @test_throws BadDataStructure project!(copy(s), P"IIZ"; keep_result=true)  # operators
        @test_throws BadDataStructure project!(copy(s), P"IZI"; keep_result=false) # in or out of
        @test_throws BadDataStructure project!(copy(s), P"IIZ"; keep_result=false) # the stabilizer
        s = S"ZII
              IZI
              IIZ"
        s = Destabilizer(s)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, 0x0)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing)
        s = S"ZII
              IZI"
        s = MixedDestabilizer(s)
        mds, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        @test mds.rank == 2
        mds, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        @test mds.rank == 3
        mds, a, r = project!(copy(s), P"IZI"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator in the stabilizer
        @test mds.rank == 2
        mds, a, r = project!(copy(s), P"IIZ"; keep_result=false)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        @test mds.rank == 3
    end
    @testset "Results from canonicalization vs from destabilizer" begin
        for n in test_sizes
            for r in [n, rand(n÷3:n*2÷3)]
                if r==0
                    @test_broken error("can not process empty stab")
                    continue
                end
                s = random_stabilizer(r,n)
                ms = MixedStabilizer(copy(s))
                d = Destabilizer(copy(s))
                md = MixedDestabilizer(copy(s))
                p = random_pauli(n,realphase=true)
                _, as, rs = project!(s,p)
                _, ams, rms = project!(ms,p)
                _, amd, rmd = project!(md,p)
                @test rs == rms == rmd
                @test canonicalize!(s) == canonicalize!(stabilizerview(ms)) == canonicalize!(stabilizerview(md))
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
        s = S"ZII
              IZI"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer

        s = S" ZII
              -IZI"
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer

        s = S"ZII
              IZI
              III"
        s = MixedStabilizer(s, 2)
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        s = S" ZII
              -IZI
               III"
        s = MixedStabilizer(s, 2)
        _, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        _, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer

        s = S"ZII
              IZI
              IIZ"
        s = Destabilizer(s)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, 0x0)
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (0, 0x2)
        s = S" ZII
               IZI
              -IIZ"
        s = Destabilizer(s)
        _, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, 0x2)
        _, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (0, 0x0)

        s = S"ZII
              IZI"
        s = MixedDestabilizer(s)
        mds, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
        mds, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        mds, a, r = project!(copy(s), P"IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        mds, a, r = project!(copy(s), P"-IIZ"; keep_result=true)
        @test (a, r) == (0, nothing) # on commuting operator out of the stabilizer
        s = S" ZII
              -IZI"
        s = MixedDestabilizer(s)
        mds, a, r = project!(copy(s), P"IZI"; keep_result=true)
        @test (a, r) == (0, 0x2) # on commuting operator in the stabilizer
        mds, a, r = project!(copy(s), P"-IZI"; keep_result=true)
        @test (a, r) == (0, 0x0) # on commuting operator in the stabilizer
    end
    @testset "Redundant row permutations in `project!(::MixedDestabilizer)`" begin
        # Fixed in 41ed1d3c
        destab =  S"+ ZX_Y_YXZ
                    + XY_Y____
                    + _Z_XXY__
                    + _ZYXXY__
                    + X__Y_ZXZ
                    + X__YXZXZ
                    + ___YXXZZ
                    + _______Z"
        stab =    S"+ X_______
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
end

if doset("Partial traces")
@testset "Partial traces" begin
    @testset "RREF canonicalization vs manual traceout" begin
        for N in test_sizes
            for n in [N,rand(N÷4:N÷2)]
                if n==0
                    @test_broken error("can not process empty stab")
                    continue
                end
                to_delete = randperm(N)[1:rand(N÷4:N÷2)]
                stab0 = random_stabilizer(n, N)
                id_paulis = zero(PauliOperator, N)
                # Trace out by doing projective measurements
                naive_stab = copy(stab0)
                for i in to_delete
                    naive_stab, anticom_index, result = project!(naive_stab, single_x(N,i))
                    if anticom_index!=0
                        naive_stab[anticom_index] = id_paulis
                    end
                    naive_stab, anticom_index, result = project!(naive_stab, single_z(N,i))
                    if anticom_index!=0
                        naive_stab[anticom_index] = id_paulis
                    end
                end
                canonicalize!(naive_stab)
                # Trace out by using the RREF canonical form
                stab = copy(stab0)
                stab, last_row = canonicalize_rref!(stab, to_delete)
                for i in last_row+1:n
                    stab[i] = id_paulis
                end
                canonicalize!(stab)
                # Confirm the results are the same
                @test stab == naive_stab
                @test mixed_stab_looks_good(stab[1:last_row])
                # Check the built-in traceout! functions for this
                s = traceout!(copy(stab0), to_delete)
                canonicalize!(s)
                @test stab == s
                # On MixedStabilizer instances
                s = traceout!(MixedStabilizer(copy(stab0), n), to_delete)
                canonicalize!(s)
                @test stab[1:last_row] == stabilizerview(s)
                @test mixed_stab_looks_good(s)
                # On MixedDestabilizer instances
                s = traceout!(MixedDestabilizer(copy(stab0)), to_delete)
                @test mixed_destab_looks_good(s)
                s = canonicalize!(stabilizerview(s))
                @test stab[1:last_row] == s
            end
        end
    end
end
end

if doset("Qubit resets")
    @testset "Qubit resets" begin
        for N in test_sizes
            for R in [rand(N÷2:N*2÷3), N]
                if N<10
                    @test_broken error("can not process empty stab")
                    continue
                end
                s = random_stabilizer(R,N)
                newstate = random_stabilizer(rand(N÷4:N*2÷3))
                perm = randperm(N)[1:nqubits(newstate)]
                to_trace = setdiff(1:N,perm)
                # Testing MixedDestabilizer
                md = MixedDestabilizer(s)
                mdr1 = reset_qubits!(copy(md), newstate,perm)
                @test mixed_destab_looks_good(mdr1)
                mdr2 = reset_qubits!(copy(mdr1),newstate,perm)
                @test mdr1==mdr2
                traceout!(mdr2,to_trace)
                mdr2v = stabilizerview(mdr2)
                @test canonicalize!(copy(mdr2v)[:,perm]) == canonicalize!(copy(newstate))
                # Testing MixedStabilizer
                ms = MixedStabilizer(s)
                msr1 = reset_qubits!(copy(ms), newstate,perm)
                @test mixed_stab_looks_good(msr1)
                msr2 = reset_qubits!(copy(msr1),newstate,perm)
                @test msr1==msr2
                traceout!(msr2,to_trace)
                msr2v = stabilizerview(msr2)
                @test canonicalize!(copy(msr2v)[:,perm]) == canonicalize!(copy(newstate))
                @test canonicalize!(msr2v) == canonicalize!(mdr2v)
                # Testing Stabilizer
                ss = R==N ? s : MixedStabilizer(s).tab # Ensure the tableau is padded with Is
                ssr1 = reset_qubits!(copy(ss), newstate,perm)
                ssr2 = reset_qubits!(copy(ssr1),newstate,perm)
                @test canonicalize!(ssr1)==canonicalize!(ssr2)
                traceout!(ssr2,to_trace)
                ssr2v = stabilizerview(ssr2)
                c, x, z = canonicalize!(ssr2v, ranks=true)
                @test canonicalize!(copy(ssr2v)[:,perm])[1:z] == canonicalize!(copy(newstate))
                @test canonicalize!(msr2v) == c[1:z]
            end
        end
    end
end    
end

function test_cliff()
if doset("Clifford Operators")
@testset "Clifford Operators" begin
    @testset "Permutations of qubits" begin
        for c in [CNOT, CliffordId⊗Hadamard, CNOT⊗CNOT, tensor_pow(CNOT,6), tensor_pow(CNOT,7), tensor_pow(CNOT,6)⊗Phase, tensor_pow(CNOT,7)⊗Phase]
            for rep in 1:5
                p = randperm(nqubits(c))
                s = random_stabilizer(nqubits(c))
                @test permute(c,p)*s[:,p] == (c*s)[:,p]
            end
        end
        for i in 1:5
            p = randperm(125)
            c = rand([CliffordId, Hadamard, Phase], 125)
            @test ⊗(c[p]...) == permute(⊗(c...), p)
        end
    end
    @testset "Tensor products" begin
        for n in test_sizes
            for np in [2,3,4]
                for pow in [1,2,10]
                    s1 = random_stabilizer(n)
                    sps = [random_stabilizer(np) for i in 1:pow]
                    ss = [s1, sps...]
                    c1 = random_clifford(n)
                    cps = repeat([random_clifford(np)],pow)
                    cs = [c1, cps...]
                    res1 = ⊗([c*s for (c,s) in zip(cs,ss)]...)
                    res2 = ⊗(cs...)*⊗(ss...)
                    res3 = (c1*s1) ⊗ (⊗(tensor_pow(cps[1],pow)) * ⊗(sps...))
                    @test res1==res2==res3
                end
            end
        end
    end
    @testset "Clifford acting on Stabilizer" begin
        for size in test_sizes
            size < 5 && continue
            s = random_stabilizer(size)
            gates = vcat([CNOT, Hadamard, Phase], repeat([CliffordId],size-4))
            gates_perm = randperm(size-1)
            gates = gates[gates_perm]
            big_gate = reduce(⊗,gates)

            s1 = apply!(copy(s),big_gate)
            @test stab_looks_good(s1)

            igates_perm = perm_inverse(gates_perm)
            s2 = copy(s)
            canonicalize!(s2)
            s2 = apply!(s2, CNOT, [igates_perm[1],igates_perm[1]+1])
            canonicalize!(s2)
            s2 = apply!(s2, Hadamard, [igates_perm[2]+(igates_perm[1]<igates_perm[2])])
            canonicalize!(s2)
            s2 = apply!(s2, Phase, [igates_perm[3]+(igates_perm[1]<igates_perm[3])])

            @test canonicalize!(s1) == canonicalize!(s2)
        end
    end
    @testset "Clifford acting on (Mixed)(De)Stabilizer" begin
        for size in test_sizes
            size < 2 && continue # TODO remove this line
            d = random_destabilizer(size)
            md = MixedDestabilizer(copy(d),size÷2)
            s = copy(stabilizerview(md))
            c = random_clifford(size)
            @test QuantumClifford._apply_nonthread!(s,c) == stabilizerview(apply!(md,c)) == stabilizerview(MixedDestabilizer(apply!(d,c),size÷2))
            newsize = min(size, 5)
            indices = randperm(size)[1:newsize]
            cn = random_clifford(newsize)
            @test QuantumClifford._apply_nonthread!(s,cn,indices) == stabilizerview(apply!(md,cn,indices)) == stabilizerview(MixedDestabilizer(apply!(d,cn,indices),size÷2))
        end
    end
    @testset "Inversions" begin
        for n in test_sizes
            c = random_clifford(n)
            ci = inv(c; phases=false)
            cip = inv(c)
            id = one(c)
            @test c*cip == cip*c == id
            @test (ci*c).tab.xzs == (c*ci).tab.xzs == id.tab.xzs
        end
    end
end
end

if doset("Clifford Operators (column representation)")
@testset "Clifford Operators (column representation)" begin
    @testset "Permutations of qubits" begin
        function naive_permute(c::CliffordColumnForm,p::AbstractArray{T,1} where T) # this is extremely slow stupid implementation
            ops = QuantumClifford.getallpaulis(c)
            CliffordColumnForm([ops[i][p] for i in 1:2*c.nqubits][vcat(p,p.+c.nqubits)])
        end
        for c in [CNOTcol, CliffordIdcol⊗Hadamardcol, CNOTcol⊗CNOTcol, tensor_pow(CNOTcol,6), tensor_pow(CNOTcol,7), tensor_pow(CNOTcol,6)⊗Phasecol, tensor_pow(CNOTcol,7)⊗Phasecol]
            for rep in 1:5
                p = randperm(nqubits(c))
               @test permute(c,p) == naive_permute(c,p)
            end
        end
        for i in 1:5
            p = randperm(125)
            c = rand([CliffordIdcol, Hadamardcol, Phasecol], 125)
            @test ⊗(c[p]...) == naive_permute(⊗(c...), p) == permute(⊗(c...), p)
        end
    end
    @testset "Tensor products" begin
        function naive_mul(l::CliffordColumnForm, r::CliffordColumnForm) # this is extremely slow stupid implementation
            opsl = QuantumClifford.getallpaulis(l)
            opsr = QuantumClifford.getallpaulis(r)
            onel = zero(opsl[1])
            oner = zero(opsr[1])
            opsl = [l⊗oner for l in opsl]
            opsr = [onel⊗r for r in opsr]
            CliffordColumnForm(vcat(opsl[1:end÷2],opsr[1:end÷2],opsl[end÷2+1:end],opsr[end÷2+1:end]))
        end
        function naive_tensor_pow(op::CliffordColumnForm,power::Integer,mem::Dict{Integer,CliffordColumnForm})
            if power==1
                return op
            elseif haskey(mem,power)
                return mem[power]
            end
            half,rest = divrem(power,2)
            phalf = get!(mem,half) do
                naive_tensor_pow(op,half,mem)
            end
            res = naive_mul(phalf,phalf)
            if rest!=0
                prest = get!(mem,rest) do
                    naive_tensor_pow(op,half,mem)
                end
                res = naive_mul(res,prest)
            end
            res
        end
        function naive_tensor_pow(op::CliffordColumnForm,power::Integer)
            naive_tensor_pow(op,power,Dict{Integer,CliffordColumnForm}())
        end
        for s in test_sizes
            for g in [CNOTcol,Hadamardcol,CNOTcol⊗Phasecol]
                @test naive_tensor_pow(g,s)==tensor_pow(g,s)
            end
        end
        for g1 in [CNOTcol,Hadamardcol,CNOTcol⊗Phasecol,CliffordIdcol]
            for g2 in [CNOTcol,Hadamardcol,CNOTcol⊗Phasecol,CliffordIdcol]
                @test naive_mul(g1,g2)==g1⊗g2
            end
        end
        c1 = tensor_pow(CNOTcol,32)
        @test naive_mul(c1,c1) == c1⊗c1
        c1 = naive_tensor_pow(Hadamardcol,33)
        @test naive_mul(c1,c1) == c1⊗c1
        c2 = naive_tensor_pow(Hadamardcol,32)
        @test naive_mul(c1,c2) == c1⊗c2
        @test naive_mul(c2,c1) == c2⊗c1
    end
    @testset "Clifford acting on Stabilizer" begin
        for size in test_sizes
            size < 5 && continue
            s = random_stabilizer(size)
            gates = vcat([CNOTcol, Hadamardcol, Phasecol], repeat([CliffordIdcol],size-4))
            gates_perm = randperm(size-1)
            gates = gates[gates_perm]
            big_gate = reduce(⊗,gates)

            s1 = apply!(copy(s),big_gate; phases=false)
            @test stab_looks_good(s1)

            igates_perm = perm_inverse(gates_perm)
            s2 = copy(s)
            s2 = apply!(s2, CNOTcol, [igates_perm[1],igates_perm[1]+1]; phases=false)
            s2 = apply!(s2, Hadamardcol, [igates_perm[2]+(igates_perm[1]<igates_perm[2])]; phases=false)
            s2 = apply!(s2, Phasecol, [igates_perm[3]+(igates_perm[1]<igates_perm[3])]; phases=false)

            @test ==(s1, s2, phases=false)
        end
    end
end
end
end

function test_symbolic()
if doset("Small symbolic operators")
    @testset "Small symbolic operators" begin
        for n in test_sizes
            for i in 1:6
                op = enumerate_single_qubit_gates(i, qubit=n, phases=rand(Bool,2))
                op0 = enumerate_single_qubit_gates(i, qubit=n) 
                op_cc = CliffordOperator(op, 1, compact=true)
                op_c = CliffordOperator(op, n)
                @test SingleQubitOperator(op)==SingleQubitOperator(op_cc, n)
                op0_c = CliffordOperator(op0, n)
                s = random_stabilizer(n)
                @test apply!(copy(s),op)==apply!(copy(s),SingleQubitOperator(op))==apply!(copy(s),op_cc,[n])==apply!(copy(s),op_c)
                @test ==(apply!(copy(s),op,phases=false),apply!(copy(s),op_cc,[n],phases=false), phases=false)
                @test apply!(copy(s),op0)==apply!(copy(s),op0_c)
            end
            i = n÷2+1
            @test apply!(copy(s),sX(i)) == apply_single_x!(copy(s),i)
            @test apply!(copy(s),sY(i)) == apply_single_y!(copy(s),i)
            @test apply!(copy(s),sZ(i)) == apply_single_z!(copy(s),i)
            n==1 && continue
            s = random_stabilizer(n)
            i1,i2 = randperm(n)[1:2]
            @test apply!(copy(s),CNOT,[i1,i2]) == apply!(copy(s),sCNOT(i1,i2))
            @test apply!(copy(s),SWAP,[i1,i2]) == apply!(copy(s),sSWAP(i1,i2))
        end
    end
    @test_throws DimensionMismatch SingleQubitOperator(CNOT,1)
    @test_throws DimensionMismatch CliffordOperator(sHadamard(5),2)
    @test_throws ArgumentError CliffordOperator(sHadamard(5),6,compact=true)
end
end

function test_random()
if doset("Random sampling of operators")
    @testset "Random sampling of operators" begin
        for n in [1, test_sizes..., 200,500]
            p = random_pauli(n)
            s = random_stabilizer(n)
            ss = random_stabilizer(rand(1:n),n)
            ms = MixedDestabilizer(ss)
            d = random_destabilizer(n)
            c = random_clifford(n)
            sq = random_singlequbitop(n)
            @test stab_looks_good(s)
            @test stab_looks_good(ss)
            @test destab_looks_good(d)
            @test mixed_destab_looks_good(ms)
            @test stab_looks_good(c*s)
            @test stab_looks_good(c*ss)
            @test destab_looks_good(c*d)
            @test mixed_destab_looks_good(c*ms)
            @test stab_looks_good(p*s)
            @test stab_looks_good(p*ss)
            @test destab_looks_good(p*d)
            @test mixed_destab_looks_good(p*ms)
            @test stab_looks_good(apply!(s,sq,phases=false))
            @test stab_looks_good(apply!(ss,sq,phases=false))
            @test destab_looks_good(apply!(d,sq,phases=false))
            @test mixed_destab_looks_good(apply!(ms,sq,phases=false))
        end
    end
end
end

function test_bitpack()
if doset("Alternative bit packing")
@testset "Alternative bit packing" begin
    for n in [1,3,5]
        N = 64*n-2
        s64 = random_stabilizer(N,N);
        phases = s64.phases;
        xzs64_colmajor = s64.xzs;
        xzs64_rowmajor = collect(s64.xzs')';
        p64 = random_pauli(N);
        c64_stab = Destabilizer(random_stabilizer(N,N)).tab
        
        after_p = stab_to_gf2(p64*s64)
        after_p_phases = (p64*s64).phases
        after_can = stab_to_gf2(canonicalize!(copy(s64)))
        after_cliff = stab_to_gf2(apply!(copy(s64),CliffordOperator(c64_stab)))

        for int in [UInt8, UInt16, UInt32, UInt64], order in [:column,:row]
            p = PauliOperator(p64.phase, N, collect(reinterpret(int,p64.xz)));
            xzs_rowmajor = collect(reinterpret(int, collect(xzs64_colmajor')))';
            xzs_colmajor = collect(xzs_rowmajor);
            s_col = Stabilizer(phases,N,xzs_colmajor);
            s_row = Stabilizer(phases,N,xzs_rowmajor);
            s = order == :column ? s_col : s_row
            apply_pauli = p*s
            @test after_p_phases == apply_pauli.phases
            canon = canonicalize!(deepcopy(s))
            @test after_can == stab_to_gf2(canon)

            for c_order in [:column, :row]
                c_raw = Stabilizer(zeros(UInt8, 2N), N, collect(reinterpret(int, collect(c64_stab.xzs'))'))
                c = CliffordOperator(c_raw)
                if c_order == :row
                    c = CliffordOperator(Stabilizer(c.tab.phases, N, collect(c.tab.xzs')'))
                end
                @test after_cliff == stab_to_gf2(apply!(deepcopy(s),c))
                c_raw = CliffordColumnForm(c_raw)
                c = CliffordColumnForm(c_raw.phases, N,
                    collect(reinterpret(int, c_raw.xztox')'),
                    collect(reinterpret(int, c_raw.xztoz')'))
                if c_order == :row
                    c = CliffordColumnForm(c.phases, N, collect(c.xztox')', collect(c.xztoz')')
                end
                @test after_cliff == stab_to_gf2(apply!(deepcopy(s),c;phases=false))
            end
        end
    end
end
end
end

function test_noisycircuits()
if doset("Noisy Circuits")
@testset "Noisy Circuits" begin

if doset("Monte Carlo sims")
@testset "Monte Carlo sims" begin
    @testset "Purification examples" begin
        g1 = SparseGate(CNOT, [1,3])
        g2 = SparseGate(CNOT, [2,4])
        m = BellMeasurement([X,X],[3,4])
        good_bell_state = S"XX
                            ZZ"
        canonicalize_rref!(good_bell_state)
        v = VerifyOp(good_bell_state,[1,2])
        n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.01))
        with_purification = mctrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v], trajectories=500)
        @test with_purification[:detected_failure] > 5
        @test with_purification[:undetected_failure] > 10
        @test with_purification[:true_success] > 430
        without_purification = mctrajectories(good_bell_state⊗good_bell_state, [n,v], trajectories=500)
        @test without_purification[:detected_failure] == 0
        @test without_purification[:undetected_failure] > 10
        @test without_purification[:true_success] > 450
        nonoise = mctrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v], trajectories=10)
        @test nonoise[:detected_failure] == 0
        @test nonoise[:undetected_failure] == 0
        @test nonoise[:true_success] == 10
    end
end
end

if doset("Perturbative expansion sims")
@testset "Perturbative expansion sims" begin
    @testset "Purification examples comparison to MC" begin
        compare(a,b, symbol) = abs(a[symbol]/500-b[symbol]) / (a[symbol]/500+b[symbol]+1e-5) < 0.3
        g1 = SparseGate(CNOT, [1,3])
        g2 = SparseGate(CNOT, [2,4])
        m = BellMeasurement([X,X],[3,4])
        good_bell_state = S"XX
                            ZZ"
        canonicalize_rref!(good_bell_state)
        v = VerifyOp(good_bell_state,[1,2])
        n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.01))
        mc = mctrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v], trajectories=500)
        pe = petrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v])
        @test compare(mc,pe,:detected_failure)
        @test compare(mc,pe,:undetected_failure)
        @test compare(mc,pe,:true_success)
        mc = mctrajectories(good_bell_state⊗good_bell_state, [n,v], trajectories=500)
        pe = petrajectories(good_bell_state⊗good_bell_state, [n,v])
        @test compare(mc,pe,:detected_failure)
        @test compare(mc,pe,:undetected_failure)
        @test compare(mc,pe,:true_success)
        mc = mctrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v], trajectories=500)
        pe = petrajectories(good_bell_state⊗good_bell_state, [g1,g2,m,v])
        @test compare(mc,pe,:detected_failure)
        @test compare(mc,pe,:undetected_failure)
        @test compare(mc,pe,:true_success)
    end
    
    @testset "Symbolic" begin
        for statetype in [Stabilizer, MixedDestabilizer]
            R, (e,) = AbstractAlgebra.PolynomialRing(AbstractAlgebra.RealField, ["e"])
            unity = R(1);
            
            good_bell_state = statetype(S"XX
                                          ZZ")
            initial_state = good_bell_state⊗good_bell_state

            g1 = SparseGate(CNOT, [1,3]) # CNOT between qubit 1 and qubit 3 (both with Alice)
            g2 = SparseGate(CNOT, [2,4]) # CNOT between qubit 2 and qubit 4 (both with Bob)
            m = BellMeasurement([X,X],[3,4]) # Bell measurement on qubit 3 and 4
            v = VerifyOp(good_bell_state,[1,2]) # Verify that qubit 1 and 2 indeed form a good Bell pair
            epsilon = e # The error rate
            n = NoiseOpAll(UnbiasedUncorrelatedNoise(epsilon))

            # This circuit performs a depolarization at rate `epsilon` to all qubits,
            # then bilater CNOT operations
            # then a Bell measurement
            # followed by checking whether the final result indeed corresponds to the correct Bell pair.
            circuit = [n,g1,g2,m,v]

            pe_symbolic = petrajectories(initial_state, circuit, branch_weight=unity) # perturbative expansion
            @test pe_symbolic[:undetected_failure] == -162.0*e^4 + 162.0*e^3 + -54.0*e^2 + 6.0*e
            @test pe_symbolic[:detected_failure]   == -108.0*e^4 + 108.0*e^3 + -36.0*e^2 + 4.0*e
            @test pe_symbolic[:true_success]       == 27.0*e^4 + -54.0*e^3 + 36.0*e^2 + -10.0*e + 1.0
        end
    end
end
end

if doset("Measurements")
@testset "Measurements" begin
    # compare(a, b, symbols, c) = [abs(a[symbol]/c-b[symbol]) / (a[symbol]/c+b[symbol]+1e-5) < 0.1 for symbol in symbols]
    @testset "BellMeasurements" begin
        stateX = S"X"
        mX = BellMeasurement([X], [1])
        vX = VerifyOp(S"X", [1])
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            determinate1 = mctrajectories(stabType(stateX), [mX,vX], trajectories=10)
            @test determinate1[:detected_failure] == 0
            @test determinate1[:undetected_failure] == 0
            @test determinate1[:true_success] == 10
            determinate1_pe = petrajectories(stabType(copy(stateX)), [mX,vX])
            @test determinate1_pe[:detected_failure] == 0
            @test determinate1_pe[:undetected_failure] == 0
            @test determinate1_pe[:true_success] == 1
        end
        stateZ = S"Z"
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            random1 = mctrajectories(stabType(stateZ), [mX,vX], trajectories=500)
            @test random1[:detected_failure] > 220
            @test random1[:undetected_failure] == 0
            @test random1[:true_success] > 220
            random1_pe = petrajectories(stabType(copy(stateZ)), [mX,vX])
            @test random1_pe[:detected_failure] > 0.44
            @test random1_pe[:undetected_failure] == 0
            @test random1_pe[:true_success] > 0.44
        end
        bell_state = S" XX
                        ZZ"
        m1 = BellMeasurement([X,X], [1,2])
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            determinate2 = mctrajectories(stabType(bell_state), [m1], trajectories=10)
            @test determinate2[:detected_failure] == 0
            @test determinate2[:undetected_failure] == 0
            @test determinate2[:continue] == 10
            determinate2_pe = petrajectories(stabType(copy(bell_state)), [m1])
            @test determinate2_pe[:detected_failure] == 0
            @test determinate2_pe[:undetected_failure] == 0
        end
        m2 = BellMeasurement([X,Z], [1,2])
        v = VerifyOp(bell_state, [1,2])
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            random2 = mctrajectories(stabType(bell_state), [m2,v], trajectories=500)
            @test random2[:detected_failure]+random2[:undetected_failure] == 500
            @test random2[:true_success] == 0
            random2_pe = petrajectories(stabType(copy(bell_state)), [m2,v])
            @test random2_pe[:detected_failure]+random2_pe[:undetected_failure] == 1
            @test random2_pe[:true_success] == 0
        end
    end
    @testset "DenseMeasurements" begin
        ghzState = S"XXX
                    ZZI
                    IZZ"
        m1 = DenseMeasurement(P"ZZI", 1)
        v = VerifyOp(ghzState, [1,2,3])
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            register1 = Register(stabType(ghzState), zeros(Bool, 1))
            determinate1 = mctrajectories(register1, [m1,v], trajectories=10)
            @test determinate1[:detected_failure] == 0
            @test determinate1[:undetected_failure] == 0
            @test determinate1[:true_success] == 10
            determinate1_pe = petrajectories(register1, [m1,v])
            @test determinate1_pe[:detected_failure] == 0
            @test determinate1_pe[:undetected_failure] == 0
            @test determinate1_pe[:true_success] == 1
        end
        m2 = DenseMeasurement(P"ZII", 1)
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            register1 = Register(stabType(ghzState), zeros(Bool, 1))
            random1 = mctrajectories(register1, [m2,v], trajectories=50)
            @test random1[:detected_failure] == 0
            @test random1[:undetected_failure] == 50
            @test random1[:true_success] == 0
            random1_pe = petrajectories(register1, [m2,v])
            @test random1_pe[:detected_failure] == 0
            @test random1_pe[:undetected_failure] == 1
            @test random1_pe[:true_success] == 0
        end
        m3 = DenseMeasurement(P"XII", 1)
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            register1 = Register(stabType(ghzState), zeros(Bool, 1))
            random2 = mctrajectories(register1, [m3,v], trajectories=50)
            @test random2[:detected_failure] == 0
            @test random2[:undetected_failure] == 50
            @test random2[:true_success] == 0
            random2_pe = petrajectories(register1, [m3,v])
            @test random2_pe[:detected_failure] == 0
            @test random2_pe[:undetected_failure] == 1
            @test random2_pe[:true_success] == 0
        end
    end
    @testset "SparseMeasurements" begin
        ghzState = S"XXX
                    ZZI
                    IZZ"
        m1 = SparseMeasurement(P"ZZ", [1,2], 1)
        v = VerifyOp(ghzState, [1,2,3])
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            register1 = Register(stabType(ghzState), zeros(Bool, 1))
            determinate1 = mctrajectories(register1, [m1,v], trajectories=10)
            @test determinate1[:detected_failure] == 0
            @test determinate1[:undetected_failure] == 0
            @test determinate1[:true_success] == 10
            determinate1_pe = petrajectories(register1, [m1,v])
            @test determinate1_pe[:detected_failure] == 0
            @test determinate1_pe[:undetected_failure] == 0
            @test determinate1_pe[:true_success] == 1
        end
        m2 = SparseMeasurement(P"Z", [1], 1)
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            register1 = Register(stabType(ghzState), zeros(Bool, 1))
            random1 = mctrajectories(register1, [m2,v], trajectories=50)
            @test random1[:detected_failure] == 0
            @test random1[:undetected_failure] == 50
            @test random1[:true_success] == 0
            random1_pe = petrajectories(register1, [m2,v])
            @test random1_pe[:detected_failure] == 0
            @test random1_pe[:undetected_failure] == 1
            @test random1_pe[:true_success] == 0
        end
        m3 = SparseMeasurement(P"X", [1], 1)
        for stabType in [Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer]
            register1 = Register(stabType(ghzState), zeros(Bool, 1))
            random2 = mctrajectories(register1, [m3,v], trajectories=50)
            @test random2[:detected_failure] == 0
            @test random2[:undetected_failure] == 50
            @test random2[:true_success] == 0
            random2_pe = petrajectories(register1, [m3,v])
            @test random2_pe[:detected_failure] == 0
            @test random2_pe[:undetected_failure] == 1
            @test random2_pe[:true_success] == 0
        end
    end
end
end

if doset("Classical Bits")
@testset "Classical Bits" begin
    @testset "DecisionGate" begin
        X_error = CliffordOperator([P"X", P"-Z"])
        # testing single digit return value from decision function
        for s in [S"Z", S"-Z", S"X", S"-X", S"Y", S"-Y"]
            r = Register(s, [false])
            applyop!(r, DenseMeasurement(P"Z", 1))
            correctiveGate = SparseGate(X_error, [1])
            decisionFunction = syndrome -> syndrome[1] ? 1 : nothing
            applyop!(r, DecisionGate([correctiveGate], decisionFunction))
            @test r.stab == S"Z"
        end

        # testing an array return from decision function
        expectedFinalState = S"ZI
                               IZ"
        s = QuantumClifford.bell()
        r = Register(s, [false])
        applyop!(r, DenseMeasurement(P"ZI", 1))
        # applyop!(r, SparseMeasurement(P"Z", [1], 1))
        correctiveGates = [SparseGate(X_error, [1]), SparseGate(X_error, [2])]
        decisionFunction = syndrome -> syndrome[1] ? [1,2] : nothing
        applyop!(r, DecisionGate(correctiveGates, decisionFunction))
        canonicalize!(r.stab)
        @test r.stab == expectedFinalState

        s = QuantumClifford.bell((false, true)) # |01>+|10>
        r = Register(s, [false])
        # applyop!(r, DenseMeasurement(P"ZI", 1))
        applyop!(r, SparseMeasurement(P"Z", [1], 1))
        # we use the same corrective gates, with a different decision function
        decisionFunction = syndrome -> syndrome[1] ? [1] : 2 # both [1] and 1 should work
        applyop!(r, DecisionGate(correctiveGates, decisionFunction))
        canonicalize!(r.stab)
        @test r.stab == expectedFinalState
    end
    @testset "ConditionalGate" begin
        id_op = CliffordOperator([P"X", P"Z"])
        X_error = CliffordOperator([P"X", P"-Z"])

        for s in [S"Z", S"-Z", S"X", S"-X", S"Y", S"-Y"]
            r = Register(s, [false])
            applyop!(r, DenseMeasurement(P"Z", 1))
            correctiveGate = SparseGate(X_error, [1])
            identityGate = SparseGate(id_op, [1])
            applyop!(r, ConditionalGate(correctiveGate, identityGate, 1))
            @test r.stab == S"Z"
        end

        expectedFinalState = S"ZI
                               IZ"
        s = QuantumClifford.bell((false, true))
        r = Register(s, [false])
        applyop!(r, DenseMeasurement(P"ZI", 1))
        correctiveGate1 = SparseGate(X_error, [1])
        correctiveGate2 = SparseGate(X_error, [2])
        applyop!(r, ConditionalGate(correctiveGate1, correctiveGate2, 1))
        canonicalize!(r.stab)
        @test r.stab == expectedFinalState
    end
end
end

if doset("Quantikz diagrams")
@testset "Quantikz diagrams" begin
    noise = UnbiasedUncorrelatedNoise(0.1)
    @test circuit2string([
            SparseGate(CNOT, [1,4]),
            SparseGate(CNOT, [3,2]),
            SparseGate(CPHASE, [1,2]),
            SparseGate(SWAP, [2,4]),
            SparseGate(CNOT*CNOT, [1,3]),
            NoiseOp(noise,[1,3]),
            NoiseOpAll(noise),
            NoisyGate(SparseGate(CNOT*CNOT, [2,4]),noise),
            ]) == "\\begin{quantikz}[transparent, row sep={0.8cm,between origins}]\n\\qw & \\ctrl{0} & \\qw & \\ctrl{0} & \\qw & \\gate[3,label style={yshift=0.2cm},disable auto height]{\\;\\;} & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\qw & \\qw\\\\\n\\qw & \\qw & \\targ{}\\vqw{0} & \\ctrl{-1} & \\swap{0} & \\linethrough & \\qw & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\gate[3,label style={yshift=0.2cm},disable auto height]{\\;\\;} & \\qw\\\\\n\\qw & \\qw & \\ctrl{-1} & \\qw & \\qw & \\qw & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\linethrough & \\qw\\\\\n\\qw & \\targ{}\\vqw{-3} & \\qw & \\qw & \\swap{-2} & \\qw & \\qw & \\gate[1,style={starburst,starburst points=7,inner xsep=-2pt,inner ysep=-2pt,scale=0.5}]{} & \\qw & \\qw\n\\end{quantikz}"
end
end

end
end

end

Random.seed!(42)
test_paulistab()
test_bitpack()
test_cliff()
test_operations()
test_symbolic()
test_random()
test_noisycircuits()
doctests()