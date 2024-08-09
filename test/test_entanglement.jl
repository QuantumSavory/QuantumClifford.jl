@testitem "Entanglement" begin
    using Graphs

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.

    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good

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

    @testset "Entanglement calculated from clipped/rref for mixed states" begin
        for n in test_sizes
            s = random_destabilizer(rand(1:n), n)
            endpoints = rand(1:n, 2)
            subsystem_range = min(endpoints...):max(endpoints...)
            @test entanglement_entropy(copy(s), subsystem_range, Val(:clip)) == entanglement_entropy(copy(s), subsystem_range, Val(:rref))
        end
    end

    @testset "Entanglement calculated from clipped/graph/rref for pure states" begin
        for n in test_sizes
            s = random_stabilizer(n)
            endpoints = rand(1:n, 2)
            subsystem_range = min(endpoints...):max(endpoints...)
            @test entanglement_entropy(copy(s), subsystem_range, Val(:clip)) == entanglement_entropy(copy(s), subsystem_range, Val(:graph)) == entanglement_entropy(copy(s), subsystem_range, Val(:rref))
        end
    end

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
    
    @testset "Mutual Information for Clifford Circuits" begin
        for n in test_sizes
            s = random_stabilizer(n)
            endpointsA = sort(rand(1:n, 2))
            subsystem_rangeA = endpointsA[1]:endpointsA[2]
            startB = rand(subsystem_rangeA)
            endB = rand(startB:n) 
            subsystem_rangeB = startB:endB
            # The mutual information `Iⁿ(𝒶, 𝒷) = Sⁿ𝒶 + Sⁿ𝒷 - Sⁿ𝒶𝒷 for Clifford circuits is non-negative since n is 1 [li2019measurement](@cite).
            @test mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:clip)) & mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:rref)) >= 0
            @test mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:graph)) >= 0
        end
    end
end
