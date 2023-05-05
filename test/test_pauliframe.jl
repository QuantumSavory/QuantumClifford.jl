using QuantumClifford

testSizes = [0, 59, 125, 256, 10^4]
# tests pauliFrameCircuitHandler() from PauliFrame.jl
@testset "Pauli Frame Cicuit Handler" begin 
    @testset "3 qubit repetition code" begin
        for offset in testSizes
            circuit = [sCNOT(1,4+offset), QuantumClifford.PauliError(2,0.75), sCNOT(2,4+offset), sCNOT(2,5+offset), sCNOT(3,5+offset), sMZ(4+offset,1), sMZ(5+offset,2)]
            ref = [0,0]
            m, f = QuantumClifford.pauliFrameCircuitHandler(5+offset,circuit,ref,100)
            # If the x component is set on the second qubit, then the fourth and fifth qubits should also have it set
            frame_index = 1
            for frame in f
                if frame[2][1]
                    @test frame[4+offset][1] && frame[5+offset][1]
                    @test sum(m[frame_index,:]) == 2 # make sure measurements reflect this
                end
                frame_index += 1
            end
        end
    end
    @testset "GHZ Circuit" begin
        ghz_circuit = [sHadamard(1), sCNOT(1,2), sCNOT(1,3), sMZ(1,1), sMZ(2,2), sMZ(3,3)]
        ref = [0,0,0]
        m, f = QuantumClifford.pauliFrameCircuitHandler(3,ghz_circuit,ref,10^6)
        total_1s = sum(m)
        @test total_1s%3 == 0 # test that most likely all measurements were 000 or 111
        @test ((total_1s/3)/(10^6) > 0.49) && ((total_1s/3)/(10^6) < 0.51) # ratio to 0 to 1 should be about 0.5
    end
end
# same as above but uses the circuitSim() method from PauliFrame.jl instead of pauliFrameCircuitHandler()
@testset "Pauli Frame Cicuit Sim" begin 
    @testset "3 qubit repetition code" begin        
        circuit = [sX(1), sX(1), sCNOT(1,4), QuantumClifford.PauliError(2,0.75), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2)]
        ref = [0,0]
        frame = QuantumClifford.PauliFrame(100, 5, ref); QuantumClifford.initZ!(frame)
        result = QuantumClifford.circuitSim(frame, circuit); f = result.frame; m = result.measurements
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
        frame = QuantumClifford.PauliFrame(10^6, 3, ref); QuantumClifford.initZ!(frame)

        f = QuantumClifford.circuitSim(frame, ghz_circuit); m = f.measurements
        total_1s = sum(m)
        @test total_1s%3 == 0 
        @test ((total_1s/3)/(10^6) > 0.49) && ((total_1s/3)/(10^6) < 0.51) 
    end
end