@testitem "backtrajectory" begin
    @testset "circuit1" begin
        circuit = [
            sHadamard(1),
            sCNOT(1, 2),
            sMZ(1),
            sMZ(2)
        ]
        results = backtrajectory!(circuit, 2)
        println("Results: ", results)
    end
end