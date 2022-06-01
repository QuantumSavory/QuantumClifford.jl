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


import Graphs

function test_entanglement()
    @testset "Entanglement calculated by clipping and graph" begin
        num_repeat = 500
        for n in [3, 4, 5, 6]
            for i_part in 1:num_repeat
                s = random_stabilizer(n)
                endpoints = rand(1:n, 2)
                leftend = min(endpoints...)
                rightend = max(endpoints...)
                onfail(@test entanglement_cont(copy(s), (leftend, rightend))==entanglement_from_graph(s, leftend:rightend)) do
                    @debug(leftend, rightend)
                    @debug s
                    graph = Graphs.Graph(s)
                    @debug collect(Graphs.edges(graph))
                end
            end  
        end
    end
end


import LinearAlgebra

function test_graph_entanglement()
    @testset "Entanglement calculated by clipping and graph" begin
        num_repeat = 200
        for n in test_sizes
            for i_part in 1:num_repeat
                s1 = random_stabilizer(n)
                graph = Graphs.Graph(s1)
                s = Stabilizer(graph)
                endpoints = rand(1:n, 2)
                leftend = min(endpoints...)
                rightend = max(endpoints...)
                adjmat = Matrix{Bool}(Graphs.adjacency_matrix(graph))
                subsystem = leftend:rightend
                other_subsystem = filter(i->!(i in collect(subsystem)), 1:Graphs.nv(graph))
                subadjmat = adjmat[subsystem,other_subsystem]
                @test entanglement_cont(copy(s), (leftend, rightend))==LinearAlgebra.rank(subadjmat)
            end  
        end
    end
end


test_clipping()
test_entanglement()
test_graph_entanglement()
