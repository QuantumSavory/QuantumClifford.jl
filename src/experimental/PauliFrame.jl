using Random
using QuantumClifford
import QuantumClifford: apply!, mctrajectory!, _div, _mod

export PauliFrame

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
struct PauliFrame{T}
    numframes::Int
    qubits::Int
    ref::BitVector
    frame::T
    measurements::Matrix{Bool}
end

function PauliFrame(numframes, qubits, ref)
    stab = zero(Stabilizer, numframes, qubits)
    frame = PauliFrame(numframes, qubits, Bool.(ref), stab, zeros(Bool, numframes, length(ref)))
    initZ!(frame)
    return frame
end


"""
This type is used to be able to write Pauli Error Channel gates like PauliError(2, 0.75).\n
In this case, a pauli error is applied on qubit 2 with probability 0.75.\n
See also [`pauliFrameCircuitHandler`](@ref), [`circuitSim`](@ref) for examples.
"""
struct PauliError
    qubit::Int
    p::Float64
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
    T = eltype(frame.frame.tab.xzs)

    @inbounds @simd for f in 1:frame.numframes
        @simd for row in 1:size(frame.frame.tab.xzs,1)÷2
            frame.frame.tab.xzs[end÷2+row,f] = rand(T)
        end
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
    i = op.qubit
    T = eltype(frame.frame.tab.xzs)
    lowbit = T(1)
    ibig = _div(T,i-1)+1
    ismall = _mod(T,i-1)
    ismallm = lowbit<<(ismall)
    ref = frame.ref[op.bit]

    frame.measurements[:,op.bit] .= (.!iszero.(frame.frame.tab.xzs[ibig,:] .& ismallm)) .⊻ ref
    return frame
end

"""
    apply!(frame::PauliFrame, op::PauliError)

Inserts a random pauli error into all frames in the provided PauliFrame instance with probabiltiy p, on bit bit_t.
See also [`PauliError`](@ref)
"""
function apply!(frame::PauliFrame, op::PauliError)
    p = op.p
    i = op.qubit
    T = eltype(frame.frame.tab.xzs)

    lowbit = T(1)
    ibig = _div(T,i-1)+1
    ismall = _mod(T,i-1)
    ismallm = lowbit<<(ismall)

    @inline @simd for f in 1:frame.numframes
        r = rand()
        if  r < p/3 # X error
            frame.frame.tab.xzs[ibig,f] ⊻= ismallm
        elseif r < 2p/3 # Z error
            frame.frame.tab.xzs[end÷2+ibig,f] ⊻= ismallm
        elseif r < p # Y error
            frame.frame.tab.xzs[ibig,f] ⊻= ismallm
            frame.frame.tab.xzs[end÷2+ibig,f] ⊻= ismallm
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
