@testitem "Entanglement" begin
    using Graphs

    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.
    using QuantumClifford
    using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good, mutual_information, entanglement_entropy

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

    @testset "Mutual information for Clifford circuits with entanglement_entropy check" begin
        using QuantumOpticsBase
        import QuantumOpticsBase: entanglement_entropy

        for n in test_sizes
            s = random_stabilizer(n)
            endpointsA = sort(rand(1:n, 2))
            subsystem_rangeA = endpointsA[1]:endpointsA[2]
            startB = rand(subsystem_rangeA)
            endB = rand(startB:n)
            subsystem_rangeB = startB:endB
            if !isempty(intersect(subsystem_rangeA, subsystem_rangeB))
                @test_throws ArgumentError mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:clip))
                @test_throws ArgumentError mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:rref))
                @test_throws ArgumentError mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:graph))
            else
                mi_clip  = mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:clip))
                mi_rref  = mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:rref))
                mi_graph = mutual_information(copy(s), subsystem_rangeA, subsystem_rangeB, Val(:graph))
                @test mi_clip == mi_rref == mi_graph
                @test mi_clip ≥ 0
                ψ = Ket(s)
                ρ = dm(ψ)
                union_range = first(subsystem_rangeA) : last(subsystem_rangeB)
                S_A = entanglement_entropy(ρ, subsystem_rangeA, entropy_vn)
                S_B = entanglement_entropy(ρ, subsystem_rangeB, entropy_vn)
                # If A ∪ B covers the full system (1:n), set S_AB = 0 to avoid an invalid full-system trace in entanglement_entropy.
                S_AB = union_range == (1:n) ? 0.0 : entanglement_entropy(ρ, union_range, entropy_vn)
                # For a pure state: I(A:B) = [S(A) + S(B) - S(A∪B)] / 2, and convert nats → bits by dividing by log(2).
                mi_indep = (S_A + S_B - S_AB) / (2 * log(2))
                @test isapprox(mi_clip, mi_indep; atol=1e-6)
            end
        end
    end
end
