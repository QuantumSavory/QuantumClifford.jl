function test_clipped()
    @testset "Clipped gauge of stabilizer states" begin
        for n in test_sizes
            s = random_stabilizer(n)
            s_clipped = copy(s)
            canonicalize_clip!(s_clipped)
            @test logdot(s, s_clipped)==0
            @test stab_looks_good(s_clipped)
            @test canonicalize!(copy(s_clipped))==canonicalize!(copy(s))
            bg = bigram(s_clipped; clip=false)
            rows, columns = size(stabilizerview(s_clipped))
            @test all(count(==(j), bg)==2 for j in 1:columns)
        end
    end
end


function test_entanglement_from_clipped()
    @testset "Entanglement calculated from clipped" begin
        for n in test_sizes
            s = random_stabilizer(n)
            endpoints = rand(1:n, 2)
            subsystem_range = min(endpoints...):max(endpoints...)
            @test entanglement_entropy(copy(s), subsystem_range, Val(:clip))==entanglement_entropy(s, subsystem_range, Val(:graph))
        end
    end
end


function test_entanglement_from_graph()
    @testset "Entanglement calculated from graph" begin
        for n in test_sizes
            s = random_stabilizer(n)
            endpoints = rand(1:n, 2)
            subsystem_range = min(endpoints...):max(endpoints...)
            @test entanglement_entropy(copy(s), subsystem_range, Val(:graph))==entanglement_entropy(s, subsystem_range, Val(:rref))
        end
    end
end


function test_entanglement_special_cases()
    @testset "Entanglement of special cases" begin
        s = S"
            + XZZ_ZZ
            + ZX_ZZ_
            + Z_XZ_Z
            + _ZZXZZ
            + ZZ_ZXZ
            + Z_ZZZX"
        subsystem = 1:3
        @test entanglement_entropy(copy(s), subsystem, Val(:clip))==2
        @test entanglement_entropy(copy(s), subsystem, Val(:graph))==2
        @test entanglement_entropy(copy(s), subsystem, Val(:rref))==2
    end
end

test_clipped()
test_entanglement_from_clipped()
test_entanglement_from_graph()
test_entanglement_special_cases()
