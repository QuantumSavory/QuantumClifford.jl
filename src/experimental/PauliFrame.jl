using Random
using QuantumClifford

"""
    PauliFrame(numframes, numqubits, reference_measurements)

The frame field holds a tableau. This tableau is not to be viewed as a normal stabilizer tableau, although
it does conjugate the same under Clifford operations. Each row in the tableau refers to a single frame. Each column a qubit.
The values stored are X,Y,Z,_ and indicate whether the corresponding pauli operation is waiting to be applied to that qubit.

The provided reference measurements will be used to create the measurements field of this struct. They will also be 
used when performing measurement. The indices of the provided measurements must match the secondary value of future measurements
when doing sMZ(measure_this_qubit, index_of_reference_provided). See also [`apply!(frame::PauliFrame, op:sMZ)`](@ref)

# Examples
```julia-repl
julia> frame = PauliFrame(numframes, numqubits, [0,1]])

julia> apply!(frame, sMZ(5,1)) 
julia> apply!(frame, sMZ(1,2)) 

The first ['sMZ']@ref measured qubit 5 and treated 0 as its reference measurement.
The second ['sMZ']@ref measured qubit 1 and treated 1 as its reference measurement.
```
"""
struct PauliFrame{T <: AbstractStabilizer}
    numframes::Int
    qubits::Int
    ref::Vector{Bool}
    frame::T 
    measurements::Matrix{Bool}

    function PauliFrame(numframes, qubits, ref)
        new{Stabilizer}(numframes, qubits, ref, zero(Stabilizer, numframes, qubits), zeros(Bool, numframes, length(ref)))
    end

    function PauliFrame(numframes, qubits, ref, T)
        new{T}(numframes, qubits, ref, zero(T, numframes, qubits), zeros(Bool, numframes, length(ref)))
    end
end

"""
This type is used to be able to write Pauli Error Channel gates like PauliError(2, 0.75).\n
In this case, a pauli error is applied on qubit 2 with probability 0.75.\n
See also [`pauliFrameCircuitHandler`](@ref), [`circuitSim`](@ref) for examples. 
"""
struct PauliError
    qubit::Int
    p::Float16 
end

"""
    initZ!(f::PauliFrame)

Inject random Z errors over all frames and qubits for the supplied PauliFrame with probability 0.5.

Calling this after initialization is essential for simulating any non-deterministic circuit.

# Examples
```julia-repl
julia> frame = PauliFrame(numframes, qubits, ref_measurements)
julia> initZ!(frame)
```
"""
function initZ!(frame::PauliFrame)
    z_index = (2+(frame.qubits-1)÷64, 2*(1+(frame.qubits-1)÷64))
    
    @inbounds @simd for f in 1:frame.numframes 
        frame.frame.tab.xzs[z_index[1]:z_index[2],f] = rand(UInt64, 1 + (frame.qubits-1)÷64, 1)
    end
    return frame
end

""" Applies a symbolic gate from QuantumClifford to the frames of the provided PauliFrame instance."""
function apply!(f::PauliFrame, op) 
    QuantumClifford._apply!(f.frame, op; phases=Val(false))
    return f
end

"""
    apply!(frame::PauliFrame, op::sMZ)

Applies [`sMZ`](@ref) to the PauliFrame instance. 

On initialization of the PauliFrame instance, a set of reference measurements must be given. 
The sMZ gate here must be of the form sMZ(measure_this_qubit, index_of_reference_provided). This method assumes all measurements are
at the end of the circuit.
See also [`pauliFrameCircuitHandler`](@ref), [`circuitSim`](@ref) for examples, and [`PauliFrame`](@ref) for more info.
"""
function apply!(frame::PauliFrame, op::sMZ)
    bit_t = op.qubit
    x_index = (bit_t-1)÷64 +1

    # Vector that represents, for each frame, whether there was an X flip on bit_t
    x_flips = .!iszero.(frame.frame.tab.xzs[x_index,:] .& 2^((bit_t-1)%64))

    ref = frame.ref[op.bit]
    frame.measurements[:,op.bit] = x_flips  .⊻ ref
    return frame
end

"""
    apply!(frame::PauliFrame, op::PauliError)

Inserts a random pauli error into all frames in the provided PauliFrame instance with probabiltiy p, on bit bit_t. 
See also [`PauliError`](@ref)
"""
function apply!(frame::PauliFrame, op::PauliError)
    p = op.p; bit_t = op.qubit
    xz_error = [(1,1+(frame.qubits-1)÷64), (2+(frame.qubits-1)÷64, 2*(1+(frame.qubits-1)÷64))] 
    frame_error = zeros(frame.numframes)
    rand!(frame_error)
 
    # logic for more than 64 qubits
    error_bit = zeros(Int64, 1+(frame.qubits-1)÷64)
    index = 1+(bit_t-1)÷64
    error_bit[index] = 2^(bit_t - 1 - 64*(index-1)) 
    
    for f in 1:frame.numframes    
        if frame_error[f] < p
            xyz = rand([1,2,3])
            if xyz == 1 # X error
                error = xz_error[1]
                frame.frame.tab.xzs[error[1]:error[2],f] .= frame.frame.tab.xzs[error[1]:error[2],f] .⊻  error_bit
            elseif xyz == 2 # Z error 
                error = xz_error[2]
                frame.frame.tab.xzs[error[1]:error[2],f] .= frame.frame.tab.xzs[error[1]:error[2],f] .⊻  error_bit
            else # Y error
                error = xz_error[1]
                frame.frame.tab.xzs[error[1]:error[2],f] .= frame.frame.tab.xzs[error[1]:error[2],f] .⊻  error_bit
                error = xz_error[2]
                frame.frame.tab.xzs[error[1]:error[2],f] .= frame.frame.tab.xzs[error[1]:error[2],f] .⊻  error_bit
            end
        end
    end
    return frame
end

"""
    pauliFrameCircuitHandler(qubits, circuit, ref_m,  numframes=1)

Simulates an entire circuit for the user, including constructing the PauliFrame object, and calling the initZ!() funtion. 

# Inputs: 
    - Number of qubits, a circuit (refer above), a vector of reference measurements, and the number of frames desired.
    
# Output 1: 
    A matrix where each row is a frame, each column is a measurement.
           M1     M2  
     F1[ M1,F1  M2,F1]
     F1[ M1,F2  M2,F2]
    
# Output 2: 
    Another output is the the Stabilizer data structure from QuantumClifford. 
    Assuming all measurements happen at the end of the circuit, it represents the pauli frame values right before measurement starts.
# Examples
```julia-repl
julia> circuit = [sX(1), sX(1), sCNOT(1,4), QuantumClifford.PauliError(2,0.75), sCNOT(2,4), sCNOT(2,5), sCNOT(3,5), sMZ(4,1), sMZ(5,2)]
julia> ref = [0,0]
julia> m, f = QuantumClifford.pauliFrameCircuitHandler(5,circuit,ref,10)
```
"""
function pauliFrameCircuitHandler(qubits, circuit, ref_m,  numframes=1)
    frame = PauliFrame(numframes, qubits, ref_m); initZ!(frame)
    for op in circuit
        apply!(frame, op)      
    end
    return  frame.measurements, frame.frame
end
""" 
mctrajectory!(state::PauliFrame, circuit)

Alternative to  [`pauliFrameCircuitHandler`](@ref). The difference is that this takes a PauliFrame object and then
    returns the modified instance after processing the ciruict. In this function, the user must first create the PauliFrame,
    and call initZ!() if they want the random Z errors on initialization. Calling initZ!() is recommended in all "real" applications.

# Examples
```julia-repl
julia> ghz_circuit = [sHadamard(1), sCNOT(1,2), sCNOT(1,3), sMZ(1,1), sMZ(2,2), sMZ(3,3)]
julia> ref = [0,0,0]
julia> frame = QuantumClifford.PauliFrame(10^6, 3, ref); QuantumClifford.initZ!(frame)
julia> f = QuantumClifford.circuitSim(frame, ghz_circuit); m = f.measurements; frame = f.frame 
```
"""
function mctrajectory!(state::PauliFrame, circuit)
    for op in circuit
        apply!(state, op)      
    end
    return state
end
