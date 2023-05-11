using QuantumClifford

# same as above but uses the mctrajectory!() method from PauliFrame.jl instead of pauliFrameCircuitHandler()
@testset "Pauli Frame Cicuit Sim" begin
    @testset "3 qubit repetition code" begin
        circuit = [sX(1), sX(1), sCNOT(1,4), PauliError(2,0.75), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2)]
        ref = [0,0]
        frame = PauliFrame(100, 5, ref)
        result = mctrajectory!(frame, circuit)
        f = result.frame; m = result.measurements
        frame_index = 1
        for frame in f
            if frame[2][1]
                @test frame[4][1] && frame[5][1]
                @test sum(m[frame_index,:]) == 2
            end
            frame_index += 1
        end
    end

    @testset "GHZ Circuit" begin
        ghz_circuit = [sHadamard(1), sCNOT(1,2), sCNOT(1,3), sMZ(1,1), sMZ(2,2), sMZ(3,3)]
        ref = [0,0,0]
        frame = PauliFrame(10^6, 3, ref)

        f = mctrajectory!(frame, ghz_circuit)
        m = f.measurements
        total_1s = sum(m)
        @test total_1s%3 == 0
        @test ((total_1s/3)/(10^6) > 0.49) && ((total_1s/3)/(10^6) < 0.51)
    end
end
