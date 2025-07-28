@testitem "test backtrajectory" begin
    
    function test_distribution(circuit, num_samples=1000)
        results_dict = Dict{Vector{Bool}, Int}()
        for _ in 1:num_samples
            result = backtrajectory(circuit)
            results_dict[bitview(result)] = get(results_dict, bitview(result), 0) + 1
        end
        results_dict, quantumstate(backtrajectory(circuit))
    end

    @testset "circuit 1" begin
        circuit = [sHadamard(1), sCNOT(1, 2), sMZ(1, 1), sMZ(2, 2)]
        results_dict, final_state = test_distribution(circuit)

        @test haskey(results_dict, Bool[0, 0])
        @test results_dict[Bool[0, 0]] ≈ 500
        @test haskey(results_dict, Bool[1, 1])
        @test results_dict[Bool[1, 1]] ≈ 500
        @test !haskey(results_dict, Bool[0, 1])
        @test !haskey(results_dict, Bool[1, 0])
    end
end


# random circuits 
# comapre with mctrajectory