# A Pauli Frame simulation can be created by interacting directly with the struct and its methods, or by calling the 
# pauliFrameCircuitHandler function, which will run an entire circuit for you, calling the appropiate methods
# and initalizing a PauliFrame struct.
#
# Currently, all 1 and 2 qubit gates from QuantumClifford work. However, only Z basis measurement is available.
# Another current assumption about measurement is that all measurements happen at the end of the circuit. 
#
# # # # # # # Notes about the xzs matrix manipulations # # # # # # # 
# For most of the operations, I manipulated the Stabilizer.tab.xzs matrix from the QuantumClifford.jl library.
# It's a 2 by f matrtix, where f is the number of frames
# If the first index is 1, then it refers to an X component, if 2 then Z. For Y, use both i.e. [1:2]
# The value of [(X or Z), frame_number] is equal to a binary number that represents whether the provided index for X or Z is "on"
#   Example: say we want to represent ZZ_Z_ . First convert to a binary string, treating the leftmost qubit as our 0 place
#            ZZ_Z_ -> 01011 = 11 (in decimal). So to set this on frame f, do frames.tab.xzs[2,f] = 11
# Using this, most of the manipulations of the frames were programmed, especially the injection of Pauli error channel errors.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
using Random
using QuantumClifford

struct PauliFrame
    numframes::Int
    qubits::Int
    ref::Vector{Bool}
    frame::QuantumClifford.Stabilizer # each row is one frame, each column is a qubit
    measurements::Matrix{Bool}

    PauliFrame(numframes, qubits, ref) = new(numframes, qubits, ref, zero(Stabilizer, numframes, qubits), zeros(Bool, numframes, length(ref)))
end

struct PauliError
    qubit::Int
    p::Float16 #probability of non identity error
end

# Call this after initialization to apply random Z errors with 50% probability on each qubit in each frame.
# This is essential for simulating non deterministic circuits
function initZ(f::PauliFrame)
    f.frame.tab.xzs[2,:] = rand(0:2^(f.qubits)-1,f.numframes,1)
end

# applies a symbolic gate from QuantumClifford
function apply!(f::PauliFrame, op) 
    QuantumClifford._apply!(f.frame, op; phases=Val(false))
end

# Measure in the Z basis. ref is the reference measurement and has value 0 or 1. bit_t is the qubit being measured.
# Returns a vector of length equal to the number of frames. Each index corresponds to a frame.
# The value of each index is the new measurement value (0 or 1) on the provided qubit.
# TODO Currently, this assumes that after measurement, that qubit will not be used. 
#   - to fix this, we could reset the measured qubit to I or Z, but what about entanglement?
# takes sMZ data structure
#   - sMZ.qubit is the qubit to be measured
#   - sMZ.bit is the refers to the place in the reference array provided at PauliFrame init
function apply!(frame::PauliFrame, op::sMZ)
    bit_t = op.qubit
    # Vector that represents, for each frame, whether there was an X flip on bit_t
    x_flips = .!iszero.(frame.frame.tab.xzs[1,:] .& 2^(bit_t-1))

    # apply flips to the reference measurement
    ref = frame.ref[op.bit]
    frame.measurements[:,op.bit] = x_flips  .⊻ ref
    frame
end

# Inserts a random pauli error into all frames in frame with probabiltiy p, on bit bit_t
function apply!(frame::PauliFrame, op::PauliError)
    p = op.p; bit_t = op.qubit
    # TODO Not 100% sure if XOR correctly represents multiplying paulis into the frame?
    xyz_error = [(1,1), (1,2), (2,2)]
    frame_error = zeros(frame.numframes)
    rand!(frame_error)
    for f in 1:frame.numframes    
        if frame_error[f] < p
            error = rand(xyz_error)
            frame.frame.tab.xzs[error[1]:error[2],f] .= frame.frame.tab.xzs[error[1]:error[2],f] .⊻ (2)^(bit_t-1)
        end
    end
    frame
end

# Simulates an entire circuit for the user. Here is sample input:
#   circuit = [sX(1), sX(1), sCNOT(1,4), QuantumClifford.PauliError(2,0.75), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2)]
#   ref = [0,0]
#
#   m, f = QuantumClifford.pauliFrameCircuitHandler(5,circuit,ref,10)
#
# Inputs: Number of qubits, a circuit (refer above), a vector of reference measurements, and the number of frames desired.
#
# p is the probability that a pauli error channel will produce a non identity error
#      -> Then after this, X,Y,Z are all equally likely.
#
# Output 1: a matrix where each row is a frame, each column is a measurement.
#        M1     M2  
#  F1[ M1,F1  M2,F1]
#  F1[ M1,F2  M2,F2]
#
# Output 2: Another output is the the Stabilizer data structure from QuantumClifford. Assuming all measurements happen
#          at the end of the circuit, it represents the pauli frame values right before measurement starts.
function pauliFrameCircuitHandler(qubits, circuit, ref_m,  numframes=1)
    frame = PauliFrame(numframes, qubits, ref_m); initZ(frame)
    for op in circuit
        apply!(frame, op)      
    end
    return  frame.measurements, frame.frame
end
# Alternative to the above. This allows greater future flexability with other frameworks
function circuitSim(state::PauliFrame, circuit)
    for op in circuit
        apply!(state, op)      
    end
    return state
end
