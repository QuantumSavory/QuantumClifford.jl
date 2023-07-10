using Random
using QuantumClifford
using Test

@testset "Column-fast vs row-fast operations" begin
    for n in [3, 300]
        for T in (Stabilizer, Destabilizer, MixedStabilizer, MixedDestabilizer)
            s = T(random_stabilizer(n,n))
            @test fastrow(copy(s)) == fastrow(fastcolumn(copy(s))) # TODO should not need to convert to the same layout for comparisons
            @test canonicalize!(fastrow(copy(s))) == fastrow(canonicalize!(fastcolumn(copy(s)))) # TODO should not need to convert to the same layout for comparisons
            c = random_clifford(n)
            layouts = []
            for layout in (fastrow, fastcolumn)
                ss = layout(copy(s))
                @test typeof(ss) == typeof(deepcopy(ss))
                apply!(ss, c)
                apply!(ss, sCNOT(1,n-1))
                push!(layouts, ss)
            end
            @test fastrow(copy(layouts[1])) == fastrow(copy(layouts[2]))
        end
    end
end
