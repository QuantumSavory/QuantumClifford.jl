TEST1 = false # uses pauliFrameCircuitHandler
TEST2 = true # uses circuitSim
using QuantumClifford

if(TEST1)
    # 3 qubit repetition code circuit -> Check that on frames wiht Y or X's on the 4th and 5th qubits that flips occured
    # Also it should be impossible with this error channel to have only one flip -> it's either two or nothing
    println("3 qubit repetition code circuit w/ two X gates applied on the first bit at the start")

    circuit = [sX(1), sX(1), sCNOT(1,4), QuantumClifford.PauliError(2,0.75), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2)]
    ref = [0,0]

    m, f = QuantumClifford.pauliFrameCircuitHandler(5,circuit,ref,10)
    println("Frames:"); show(IOContext(stdout::IO, :limit => true), f); print("\nMeasurements\n", m)

    # Showing the random Z errors model non-deterministic circuits.
    # Check that results are either 111 or 000. The ratio should be about 1/2
    println("\nGHZ Circuit")
    ghz_circuit = [sHadamard(1), sCNOT(1,2), sCNOT(1,3), sMZ(1,1), sMZ(2,2), sMZ(3,3)]
    ref = [0,0,0]

    m, f = QuantumClifford.pauliFrameCircuitHandler(3,ghz_circuit,ref,10^6);
    println("First 10 frames measurements out 10^6 frames: ", m[1:10,:])
    print("Ratio of 000 measurements to 111 measurements: ", (sum(m)/3)/(10^6))
end
if(TEST2)
    # 3 qubit repetition code circuit -> Check that on frames wiht Y or X's on the 4th and 5th qubits that flips occured
    # Also it should be impossible with this error channel to have only one flip -> it's either two or nothing
    println("3 qubit repetition code circuit w/ two X gates applied on the first bit at the start")
    
    circuit = [sX(1), sX(1), sCNOT(1,4), QuantumClifford.PauliError(2,0.75), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2)]
    ref = [0,0]
    frame = QuantumClifford.PauliFrame(10, 5, ref); QuantumClifford.initZ(frame)

    f = QuantumClifford.circuitSim(frame, circuit)
    println("Frames:"); show(IOContext(stdout::IO, :limit => true), f.frame); print("\nMeasurements\n", f.measurements)

    # Showing the random Z errors model non-deterministic circuits.
    # Check that results are either 111 or 000. The ratio should be about 1/2
    println("\nGHZ Circuit")
    ghz_circuit = [sHadamard(1), sCNOT(1,2), sCNOT(1,3), sMZ(1,1), sMZ(2,2), sMZ(3,3)]
    ref = [0,0,0]
    frame = QuantumClifford.PauliFrame(10^6, 3, ref); QuantumClifford.initZ(frame)

    f = QuantumClifford.circuitSim(frame, ghz_circuit);
    println("First 10 frames measurements out 10^6 frames: ", f.measurements[1:10,:])
    print("Ratio of 000 measurements to 111 measurements: ", (sum(f.measurements)/3)/(10^6))
end