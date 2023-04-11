# Runs any number of  pauli frames on a given circuit
# Sample use: 
# circuit = [("CNOT", [1,4]), ("PE", [2]),("CNOT", [2,4]),("CNOT", [2,5]),("CNOT", [3,5]), ("MZ", [4]),("MZ", [5])]
# ref = [0,0];
# pauli_frames(5,circuit,ref,5)
# TODO change this input to reflect the Tableau in the rest of Quantum Clifford
# TODO make the Pauli Error channel take a parameter to determine the probability of non identity error 
# TODO various speedups

# Inputs: Number of qubits, array of operations (gate/measure/error, bits), array of reference measurements, 
# number of frames to simulate, Boolean flag for output suppression

# Outputs: An 𝐹×𝑀 matrix, where 𝐹 is the number of frames and 𝑀 was the number of measurements

function pauli_frames(qubits, circuit,ref_m, numframes=1, showFrame=true)
    # Copying these because these variables will be modified later
    circuit = deepcopy(circuit)
    ref_m = deepcopy(ref_m)
    
    # Need number of measurments for reshaping returnable
    num_m = length(ref_m)
    
    # Macros
    X = 1; Z=2
    
    # Declaration of the frame data structure
    # It's size is the number of qubits by 2 by the number of frames
    frame = zeros(Bool, qubits, 2, numframes)
    
    # randomly apply Z error with prob 1/2 on init
    # Maybe there's a faster tensor operation to do the same thing?
    for i in 1:numframes
        for n in 1:qubits
                frame[n, Z, i]  =   rand(Bool)
        end
    end
    if showFrame
        println("\nInitial Pauli frame\n", frame)
    end
    
    # It's nice to confirm the inputs, no?
    println("Circuit\n", circuit, "\n")
    println("Reference measurements\n", ref_m, "\n")
    
    # This is the returnable
    sim_m = []
    
    # Main loop for moving the pauli frame through the circuit
    while !isempty(circuit)
        # Get the next operation (gate or measurement or error channel)
        op = popfirst!(circuit)
        println(op)
        
        # op was CNOT
        if op[1]=="CNOT"
            bit_c = op[2][1]
            bit_t = op[2][2]

            frame[bit_t,X, :] .=  frame[bit_t,X,:]  .⊻ frame[bit_c,X,:] # This  ⊻  is XOR
            
        # op was Pauli error channel
        # TODO add parameter that determines the probability of non I error
        # Currently, {I,X,Y,Z} are all equally likely
        elseif op[1]=="PE"
            bit_t = op[2]
            # There's probably a faster way to do this than a for loop
            #  Not sure if XOR correctly represents multiplying paulis into the frame? 
            for i in 1:numframes
                frame[bit_t, :, i]  =  [rand(Bool) ⊻ frame[bit_t, X, i][1]  , rand(Bool)⊻ frame[bit_t, Z, i][1]]
            end
            
        # op was measurement in Z basis
        elseif op[1]=="MZ"
            ref = popfirst!(ref_m)
            bit_t = op[2]
            # add the simulated measurement to the returnable
            append!(sim_m, frame[bit_t,X,:]  .⊻ ref[1])
            
            # Reset target bit's frame - the Z component is randomly 0 or 1
            # Potential for speedup here I think
            # TODO Commented out for now because it makes reading the final frame more confusing for debugging
            #for i in 1:numframes
            #    frame[bit_t, :, i]  =  [0  , rand(Bool)]
            #end
            
        #  Hadamard swaps X and Z
        elseif op[1]== "H"
            bit_t = op[2]
            temp = frame[bit_t,X, :]
            frame[bit_t,X, :] = frame[bit_t,Z, :]
            frame[bit_t,Z, :] = temp
        end
        
    end
    
    if showFrame
       println("\nPauli frame after at end of circuit\n", frame)
    end
 
    # return matrix: each row is a frame, each column is a measurement
    # 
    # M1F1 M2F1
    # M1F2 M2F2
    
    return reshape(sim_m, (numframes, num_m))
end


