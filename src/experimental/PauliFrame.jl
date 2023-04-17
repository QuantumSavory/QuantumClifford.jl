# Runs any number of  pauli frames on a given circuit
# Inputs: Number of qubits, a circuit (refer to examples at the bottom of this file), number of frames desired,
#         a probalitiy p, and a boolean flag that supresses some print states when turned to false
#
# p is the probability that a pauli error channel will produce a non identity error
#      -> Then after this, the three Pauli operators are all equally likely.

# # # # # # # Notes about the xzs matrix manipulations # # # # # # # 
# For most of the operations, I manipulated the Stabilizer.tab.xzs matrix from the QuantumClifford.jl library.
# It's a 2 x f matrtix, where f is the number of frames
# If the first index is 1, then it refers to an X component, if 2 then Z. For Y, use both i.e. [1:2]
# The value of [(X or Z), frame] is equal to a binary number that represents whether the provided index for X or Z is "on"
#   Example: say we want to represent ZZ_Z_ . First convert to a binary string, treating the leftmost qubit as our 0 place
#            ZZ_Z_ -> 01011 = 11 (in decimal). So to set this on frame f, do frames.tab.xzs[2,f] = 11
# Using this, most of the manipulations of the frames were programming, specifically related to injecting Pauli error channel errors
#
# By using this data structure from QuantumClifford.jl, the conjugations of the pauli errors in the frames are done by using QuantumClifford.apply!()

using Random
using QuantumClifford

runTests = true
function pauli_frames(qubits, circuit,ref_m, numframes=1, p=0.25, showFrame=true)
    # Copying these because these variables will be modified later
    circuit = deepcopy(circuit)
    ref_m = deepcopy(ref_m)
    # Need number of measurements for reshaping returnable
    num_m = length(ref_m)
    
    # Declaration of the frame data structure
    frame = zero(Stabilizer, numframes, qubits)
    # Randomly apply Z error with prob 1/2 on init
    # This generates a random vector of values from 0 to the maximum binary number for the number of qubits.
    # For example, with 2 qubits: this will be 0 through 3 or 00 01 10 11. The vector length is numframes.
    frame.tab.xzs[2,:] = rand(0:2^(qubits)-1,numframes,1)
    
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
            # Not 100% sure if XOR correctly represents multiplying paulis into the frame? 
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
       println("\nPauli frame at end of circuit\n", frame)
    end
 
    #  return matrix: each row is a frame, each column is a measurement. The measurements refer to 
    # whether the result of the corresponding ith reference measurement flipped
    #               M1        M2  
    #   F1[   M1,F1  M2,F1]
    #  F1[  M1,F2   M2,F2]
    # 
    return reshape(sim_m, (numframes, num_m)), frame
end


if runTests
    if true
        # 3 bit repetition code example - general demonstration
        println("\3 qubit rep code Circuit")
        circuit = [(:"CNOT", [1,4]), (:"PE", [2]),(:"CNOT", [2,4]),(:"CNOT", [2,5]),(:"CNOT", [3,5]), (:"MZ", [4]),(:"MZ", [5])]
        ref = [0,0];

        m, f = pauli_frames(5,circuit,ref,5, 0.1)
        print(m)
    end

    if true
        # Showing the random Z errors model non-deterministic circuits.
        println("\nGHZ Circuit")
        ghz_circuit = [(:"H",[1]), (:"CNOT",[1,2]),(:"CNOT",[1,3]), (:"MZ", [1]), (:"MZ", [2]), (:"MZ", [3])]
        ref = [0,0,0]

        m, f = pauli_frames(3,ghz_circuit,ref,10^6, 0.1,false);
        println("First 10 frames measurements out 10^6 frames: ", m[1:10,:])
        println("Ratio of 000 measurements to 111 measurements: ", (sum(m)/3)/(10^6))
    end
end



