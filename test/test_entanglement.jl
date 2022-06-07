function test_clipping()
    @testset "Clipped gauge of stabilizer states" begin
        for n in test_sizes
            num_repeat = 100
            for i_repeat in 1:num_repeat
                s = random_stabilizer(n)
                s_clipped = copy(s)
                canonicalize_clip!(s_clipped)
                @test logdot(s, s_clipped)==0
                @test stab_looks_good(s_clipped)
                bg = bigram(s_clipped; clip=false)
                rows, columns = size(stabilizerview(s_clipped))
                @test all(count(==(j), bg)==2 for j in 1:columns)
            end
        end
    end
end


import LinearAlgebra

# by Krastanov, here for testing, to be move to src file
function entanglement_entropy_traceout(state, qubits_to_be_deleted)
    nb_of_qubits = nqubits(state)
    nb_of_deletions = length(qubits_to_be_deleted)
    new_state = traceout!(MixedDestabilizer(state), qubits_to_be_deleted) # TODO: this can be written better
    rank_after_deletion = LinearAlgebra.rank(new_state)
    return nb_of_qubits - rank_after_deletion - nb_of_deletions
end

import Graphs

function test_entanglement_from_clipping()
    @testset "Entanglement calculated from clipping" begin
        num_repeat = 500
        for n in test_sizes
            for i_repeat in 1:num_repeat
                s = random_stabilizer(n)
                endpoints = rand(1:n, 2)
                subsystem_range = min(endpoints...):max(endpoints...)
                onfail(@test entanglement_entropy(copy(s), subsystem_range, Val(:clipping))==entanglement_entropy_traceout(copy(s), subsystem_range)) do
                    @debug subsystem_range 
                    @debug s
                    graph = Graphs.Graph(s)
                    @debug collect(Graphs.edges(graph))
                end
            end  
        end
    end
end


function test_entanglement_from_graph()
    @testset "Entanglement calculated from graph" begin
        num_repeat = 500
        for n in test_sizes
            for i_repeat in 1:num_repeat
                s = random_stabilizer(n)
                endpoints = rand(1:n, 2)
                subsystem_range = min(endpoints...):max(endpoints...)
                onfail(@test entanglement_entropy(copy(s), subsystem_range, Val(:graph))==entanglement_entropy_traceout(copy(s), subsystem_range)) do
                    @debug subsystem_range 
                    @debug s
                    graph = Graphs.Graph(s)
                    @debug collect(Graphs.edges(graph))
                end
            end  
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
        @test entanglement_entropy(copy(s), subsystem, Val(:clipping))==2
        @test entanglement_entropy(copy(s), subsystem, Val(:graph))==2
    end
end

test_clipping()
test_entanglement_from_clipping()
test_entanglement_from_graph()
test_entanglement_special_cases()
