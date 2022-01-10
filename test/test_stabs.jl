function test_stabs()
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
                @test mixed_destab_looks_good(md)
                canonicalize!(s)
                s1 = copy(stabilizerview(md))
                canonicalize!(s1)
                @test s1 == s
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

test_stabs()