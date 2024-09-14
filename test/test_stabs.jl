@testitem "Stabilizers" begin
    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.
    @testset "Pure and Mixed state initialization" begin

        @testset "Destabilizer initialization" begin
            for n in test_sizes
                s = random_stabilizer(n)
                d = Destabilizer(s)
                @test destab_looks_good(d)
                canonicalize!(s)
                s1 = copy(stabilizerview(d))
                canonicalize!(s1)
                @test s1 == s
            end
        end
        @testset "Mixed destabilizer initialization" begin
            for n in test_sizes
                n<10 && continue
                s = random_stabilizer(rand(n÷2+1:n-4),n)
                md = MixedDestabilizer(s)
                ms = MixedStabilizer(s)
                @test mixed_stab_looks_good(ms)
                @test mixed_destab_looks_good(md)
                canonicalize!(s)
                s1 = copy(stabilizerview(md))
                canonicalize!(s1)
                @test s1 == s == stabilizerview(canonicalize!(ms))
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
                @test MixedDestabilizer(s) == md
            end
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
            stabs = [s[1:i] for s in [random_stabilizer(n) for n in [32,16,16,64,63,65,129,128,127]] for i in rand(1:10)];
            mdstabs = MixedDestabilizer.(stabs);
            @test canonicalize!(⊗(stabs...)) == canonicalize!(stabilizerview(⊗(mdstabs...)))
        end
    end

    @testset "Stabilizer indexing" begin
        s = random_stabilizer(9,10)
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
        @test axes(s) == (axes(s,1), axes(s,2)) == (Base.OneTo(9),Base.OneTo(10))
        ms = MixedStabilizer(s)
        mds = MixedStabilizer(s)
        @test length(mds) == length(ms) == 10
        @test length(s) == 9
        for n in test_sizes
            s = random_stabilizer(n)
            r1 = rand(1:n)
            ri1 = deleteat!(collect(1:n),r1)
            s1a = QuantumClifford.remove_column!(copy(s),r1)
            s1b = copy(s)[:,ri1]
            r2 = min(n,rand([63,64,65]))
            ri2 = deleteat!(collect(1:n),r2)
            s2a = QuantumClifford.remove_column!(copy(s),r2)
            s2b = copy(s)[:,ri2]
            @test stab_to_gf2(s1a) == stab_to_gf2(s1b)
            @test stab_to_gf2(s2a) == stab_to_gf2(s2b)
        end
    end

    @testset "horizontal concatenation" begin
        @test hcat(ghz(2), ghz(2)) == S"XXXX ZZZZ"
        s1 = S"YZ -XX"
        s2 = S"-ZY -YX"
        @test hcat(copy(s1), copy(s2)) == S"-YZZY XXYX"
        @test hcat(copy(s1), copy(s2), copy(s1), copy(s2)) == S"YZZYYZZY XXYXXXYX"
        @test_throws ArgumentError hcat(copy(s1), random_stabilizer(3))
        @test hcat(copy(tab(s1)), copy(tab(s2))) == T"-YZZY XXYX"
        @test hcat(copy(tab(s1)), copy(tab(s2)), copy(tab(s1)), copy(tab(s2))) == T"YZZYYZZY XXYXXXYX"
    end
end
