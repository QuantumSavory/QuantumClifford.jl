@testitem "Decoding interfaces" begin
    using QECCore
    using Test

    @testset "FactoredBitNoiseModel" begin
        em1 = FactoredBitNoiseModel([0.4,0.5])
        @test em1.probabilities == Dict([1] => [0.6,0.4], [2] => [0.5,0.5])
        @test em1.num_bits == 2
        @test isvector(em1)
        @test isindependent(em1)

        em2 = FactoredBitNoiseModel(3, Dict([1,2] => [1.0 0.0;0.0 0.0], [3] => [0.0, 1.0]))
        @test !isvector(em2)
        @test isindependent(em2)

        em3 = FactoredBitNoiseModel(3, Dict([1,2] => [1.0 0.0;0.0 0.0], [2,3] => [1.0 0.0;0.0 0.0]))
        @test !isvector(em3)
        @test !isindependent(em3)

        ep4 = depolarization_error_model(0.1, 3)
        @test ep4.probabilities == Dict([1,4] => [0.9 0.1/3;0.1/3 0.1/3], [2,5] => [0.9 0.1/3;0.1/3 0.1/3], [3,6] => [0.9 0.1/3;0.1/3 0.1/3])
        @test ep4.num_bits == 6
        @test !isvector(ep4)
        @test isindependent(ep4)

        @test_throws AssertionError FactoredBitNoiseModel(3, Dict([1,2,3] => [0.6,0.4,0.5]))
        @test_throws AssertionError FactoredBitNoiseModel(2, Dict([1] => [0.6,0.4,0.5],[2] => [0.5,0.5]))
        @test_throws AssertionError FactoredBitNoiseModel(2, Dict([1] => [0.6,0.4],[3] => [0.5,0.5]))
        @test_throws AssertionError FactoredBitNoiseModel(2, Dict([1] => [0.6,0.2],[2] => [0.5,0.5]))
        @test_throws AssertionError FactoredBitNoiseModel(3, Dict([1] => [0.6,0.4],[3] => [0.5,0.5]))
    end

    @testset "IndependentVectorSampler" begin
        em = FactoredBitNoiseModel(fill(0.1, 100))
        ep = sample(em, 10000, IndependentVectorSampler())
        @test size(ep) == (100, 10000)
        @test count(ep)/100/10000 ≈ 0.1 atol = 0.01
    end

    @testset "DetectorModelProblem" begin
        c = Steane7()
        pm = parity_matrix(c)
        em = FactoredBitNoiseModel(fill(0.1, 14))
        logical_matrix = fill(false, 2, 14)
        logical_matrix[1,1] = true
        logical_matrix[2,1+7] = true
        logical_matrix[1,2] = true
        logical_matrix[2,2+7] = true
        logical_matrix[1,3] = true
        logical_matrix[2,3+7] = true

        problem = DetectorModelProblem(em, pm, logical_matrix)

        @test_throws AssertionError DetectorModelProblem(em, pm, fill(false, 2, 13))
        @test_throws AssertionError DetectorModelProblem(em, pm[:,1:13], logical_matrix)
        @test_throws AssertionError DetectorModelProblem(FactoredBitNoiseModel(fill(0.1, 15)), pm, logical_matrix)

        em = FactoredBitNoiseModel(fill(0.3, 14))
        ep = sample(em, 10, IndependentVectorSampler())

        syndrome = QECCore.measure_syndrome(problem, ep)
        @test size(syndrome) == (6, 10)

        ep = fill(false, 14, 1)
        ep[7] = true # Z error on qubit 7
        @test QECCore.measure_syndrome(problem, ep) == Bool[true; true; true; false; false; false;;]

        ep2 = copy(ep)
        @test !(QECCore.check_decoding_result(problem, Bool.(mod.(ep2.-ep, 2)))[])

        ep2[1] = true
        @test (QECCore.check_decoding_result(problem, Bool.(mod.(ep2.-ep, 2)))[])

        # applying a logical operator will lead to a logical error
        ep[2] = true
        ep[3] = true
        @test (QECCore.check_decoding_result(problem, Bool.(mod.(ep2.-ep, 2)))[])
        ep2 = copy(ep)
        ep2[8:10] .= true
        @test (QECCore.check_decoding_result(problem, Bool.(mod.(ep2.-ep, 2)))[])

        # applying a stabilizer will not lead to a logical error
        ep2 = copy(ep)
        ep2[4:6] .= true
        ep2[7] = false
        @test !(QECCore.check_decoding_result(problem, Bool.(mod.(ep2.-ep, 2)))[])

        ep2 = copy(ep)
        ep2[11:14] .= true
        @test !(QECCore.check_decoding_result(problem, Bool.(mod.(ep2.-ep, 2)))[])

        samples = sample(problem, 1000, IndependentVectorSampler())
        @test decoding_error_rate(problem, samples, MatrixDecodingResult(samples.physical_bits)) == 0.0
    end
end
