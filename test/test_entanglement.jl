function test_clipping()
    @testset "Clipped gauge of stabilizer states" begin
        for n in test_sizes
            num_repeat = 100
            for i_repeat in 1:num_repeat
                s = random_stabilizer(n)
                s_clipped = copy(s)
                clip!(s_clipped)
                @test logdot(s, s_clipped)==0
                bigram = get_bigram(s_clipped; do_clip=false)
                rows, columns = size(stabilizerview(s_clipped))
                @test all(count(==(j), bigram)==2 for j in 1:columns)
            end
        end
    end
end


# functions to show info when fails
onfail(body, _::Test.Pass) = nothing
onfail(body, _::Test.Fail) = body()
onfail(body, _::Tuple{Test.Fail,T}) where {T} = body()

import LinearAlgebra

# by Krastanov, here for testing
function entropy(state, qubits_to_be_deleted)
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
            for i_part in 1:num_repeat
                s = random_stabilizer(n)
                endpoints = rand(1:n, 2)
                leftend = min(endpoints...)
                rightend = max(endpoints...)
                entanglement_new = entropy(s, leftend:rightend)
                onfail(@test entanglement_cont(copy(s), (leftend, rightend))==entanglement_new) do
                    @debug(leftend, rightend)
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
        for n in [3, 4, 5, 6]
            for i_part in 1:num_repeat
                s = random_stabilizer(n)
                endpoints = rand(1:n, 2)
                leftend = min(endpoints...)
                rightend = max(endpoints...)
                entanglement_new = entropy(s, leftend:rightend)
                onfail(@test entanglement_from_graph(s, leftend:rightend)==entanglement_new) do
                    @debug(leftend, rightend)
                    @debug s
                    graph = Graphs.Graph(s)
                    @debug collect(Graphs.edges(graph))
                end
            end  
        end
    end
end

import HDF5

function test_with_qiskit()
    @testset "Compare entanglement with qiskit data" begin
        data_path = "qiskit_data/stab_ent.h5"
        test_sizes_qiskit = [3, 4, 6, 8]
        num_repeat = 200
        HDF5.h5open(data_path, "r") do f
            for test_size in test_sizes_qiskit
                grp = f["size$test_size"]
                    for i_repeat in 1:num_repeat
                        # load from qiskit data
                        xzs = convert(Matrix{Bool}, grp["xzs"][:,:,i_repeat]) # Bool is not supported in HDF5.jl
                        xzs = transpose(xzs)
                        phases = grp["phases"][:,i_repeat]
                        leftend = grp["leftend"][i_repeat] + 1
                        rightend = grp["rightend"][i_repeat] + 1
                        entanglement_qiskit = grp["entanglement"][i_repeat]
                        # calculate by our approaches
                        s = Stabilizer(phases, xzs)
                        entanglement_clipped = entanglement_cont(copy(s), (leftend, rightend))
                        entanglement_graph = entanglement_from_graph(s, leftend:rightend)
                        onfail(
                            @test isapprox(entanglement_qiskit, entanglement_clipped, atol=1e-3) &&
                            isapprox(entanglement_qiskit, entanglement_graph, atol=1e-3)
                        ) do
                            @debug entanglement_qiskit
                            @debug entanglement_clipped
                            @debug entanglement_graph
                            @debug(leftend, rightend)
                            @debug s
                            graph = Graphs.Graph(s)
                            @debug collect(Graphs.edges(graph))
                        end
                        
                    end
            end
        end
    end
end

test_clipping()
test_entanglement_from_clipping()
test_entanglement_from_graph() #TODO: rarely failed

# uncomment the line below for cross checking with Qiskit, first generate data by https://github.com/royess/test-stab-entanglement/blob/main/get_data.py
# test_with_qiskit() 
