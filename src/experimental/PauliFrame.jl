# Runs any number of  pauli frames on a given circuit
# p is the probability that a pauli error channel will produce a non identity error
#      -> Then after this, the three Pauli operators are all equally likely.
using Random
using QuantumClifford
function pauli_frames(qubits, circuit,ref_m, numframes=1, p=0.25, showFrame=true)
    # Copying these because these variables will be modified later
    circuit = deepcopy(circuit)
    ref_m = deepcopy(ref_m)
    # Need number of measurements for reshaping returnable
    num_m = length(ref_m)
    
    # Declaration of the frame data structure
    frame = zero(Stabilizer, numframes, qubits)
    # Randomly apply Z error with prob 1/2 on init
    frame.tab.xzs[2,:] = rand(0:2^(qubits)-1,qubits,1)
    
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
        if op[1]==:"CNOT"
            bit_c = op[2][1]
            bit_t = op[2][2]
            apply!(frame, tCNOT,  [bit_c, bit_t])
            
        # op was Pauli error channel
        elseif op[1]==:"PE"
            bit_t = op[2][1]
            # There's probably a faster way to do this than a for loop?
            # The difficulty is that I need to randomly choose indices for each frame
            # Not sure if XOR correctly represents multiplying paulis into the frame? 
            xyz_error = [(1,1), (1,2), (2,2)]
            frame_error = zeros(numframes)
            rand!(frame_error)
            for f in 1:numframes    
                if frame_error[f] < p
                    error = rand(xyz_error)
                    frame.tab.xzs[error[1]:error[2],f] .= frame.tab.xzs[error[1]:error[2],f] .⊻ (2)^(bit_t-1)
                end
            end
            
        # op was measurement in Z basis
        elseif op[1]==:"MZ"
            ref = popfirst!(ref_m)
            bit_t = op[2][1]
            
            # Vector that represents, for each frame, whether there was an X flip on bit_t
            x_flips = .!iszero.(frame.tab.xzs[1,:] .& 2^(bit_t-1))
            
            # add the simulated measurement to the returnable
            append!(sim_m, x_flips  .⊻ ref[1])
            
            # TODO
            # Reset target bit's frame - the Z component is randomly 0 or 1
            # For now, we assume nothing else will happen after a qubit is measured
            
        #  Hadamard 
        elseif op[1]== :"H"
            bit_t = op[2]
            apply!(frame, tHadamard,  bit_t)
        end
        
    end
    
    if showFrame
       println("\nPauli frame after at end of circuit\n", frame)
    end
 
    #  return matrix: each row is a frame, each column is a measurement. The measurements refer to 
    # whether the result of the corresponding ith reference measurement flipped
    #               M1        M2  
    #   F1[   M1,F1  M2,F1]
    #  F1[  M1,F2   M2,F2]
    # 
    return reshape(sim_m, (numframes, num_m)), frame
end

circuit = [(:"CNOT", [1,4]), (:"PE", [2]),(:"CNOT", [2,4]),(:"CNOT", [2,5]),(:"CNOT", [3,5]), (:"MZ", [4]),(:"MZ", [5])]
ref = [0,0];

m, f = pauli_frames(5,circuit,ref,5, 0.1)
print(m)