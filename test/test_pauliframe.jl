TEST1 = true
PE = :PE
MZ = :MZ

using QuantumClifford

if(TEST1)
    # 3 qubit repetition code circuit -> Check that on frames wiht Y or X's on the 4th and 5th qubits that flips occured
    # Also it should be impossible with this error channel to have only one flip -> it's either two or nothing
    println("3 qubit rep code Circuit")

    # TODO change things to this representation circuit = [sCNOT(1,2), sX(3), MX(1),...]
    circuit = [(sCNOT, [1,4]), (PE, [2]),(sCNOT, [2,4]),(sCNOT, [2,5]),(sCNOT, [3,5]), (MZ, [4]),(MZ, [5])]
    #circuit = [(sCNOT, [1,4]), (PE, [2]), (MZ, [4])]
    ref = [0, 0]

    m, f = QuantumClifford.pauliFrameCircuitHandler(5,circuit,ref,10)
    println("Frames:"); show(IOContext(stdout::IO, :limit => true), f); print("\nMeasurements\n", m)

    # Showing the random Z errors model non-deterministic circuits.
    # Check that results are either 111 or 000 
    # The ratio should be about 1/2
    println("\nGHZ Circuit")
    ghz_circuit = [(sHadamard,[1]), (sCNOT,[1,2]),(sCNOT,[1,3]), (MZ, [1]), (MZ, [2]), (MZ, [3])]
    ref = [0,0,0]

    m, f = QuantumClifford.pauliFrameCircuitHandler(3,ghz_circuit,ref,10^6);
    println("First 10 frames measurements out 10^6 frames: ", m[1:10,:])
    println("Ratio of 000 measurements to 111 measurements: ", (sum(m)/3)/(10^6))

    # Using s instead of t gates  on 3 qubit rep code circuit
    println("\nTesting sGate on 3 qubit rep code Circuit")
    s_gate_circuit = [(sX,[1]),(sX,[1]),(sCNOT, [1,4]), (PE, [2]),(sCNOT, [2,4]),(sCNOT, [2,5]),(sCNOT, [3,5]), (MZ, [4]),(MZ, [5])]
    ref = [0,0]
    
    m, f = QuantumClifford.pauliFrameCircuitHandler(5,s_gate_circuit,ref,10)
    println("Frames:"); show(IOContext(stdout::IO, :limit => true), f); print("\nMeasurements\n", m)

end