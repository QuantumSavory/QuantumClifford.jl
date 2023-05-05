using QuantumClifford
using BenchmarkTools

const SUITE = BenchmarkGroup()
const NUMFRAMES = 1

# Generate testing circuit of an initial state corresponding the stabilizer with all X's on the diagonal and indenities elsewhere
# Then,  a CNOT is performed between each qubit and a single ancillary qubit.
# The hadamard gates are to create the initial state.
function x_diag_circuit(CIRCUIT_SIZE)
    circuit = Vector{QuantumClifford.AbstractOperation}()
    for i in 1:CIRCUIT_SIZE
        push!(circuit, sHadamard(i))
        gate = sCNOT(i, CIRCUIT_SIZE+1)
        push!(circuit, gate)
    end
    push!(circuit, sMZ(CIRCUIT_SIZE+1,1))
    return circuit
end

# Run pauliframe simulation
const ref_m = [0] 

SUITE["pauliframe"] = BenchmarkGroup(["pauliframe"])


SUITE["pauliframe"]["pauliframesim"] = BenchmarkGroup(["pauliframesim"])
SUITE["pauliframe"]["pauliframesim"]["100"] = @benchmarkable QuantumClifford.pauliFrameCircuitHandler(101,circuit,ref_m,NUMFRAMES) setup=(circuit = copy(x_diag_circuit(100)))
SUITE["pauliframe"]["pauliframesim"]["1000"] = @benchmarkable QuantumClifford.pauliFrameCircuitHandler(1001,circuit,ref_m,NUMFRAMES) setup=(circuit = copy(x_diag_circuit(1000)))
SUITE["pauliframe"]["pauliframesim"]["1000_1000"] = @benchmarkable QuantumClifford.pauliFrameCircuitHandler(1001,circuit,ref_m,1000) setup=(circuit = copy(x_diag_circuit(1000)))
SUITE["pauliframe"]["pauliframesim"]["10^6_100"] = @benchmarkable QuantumClifford.pauliFrameCircuitHandler(10^6+1,circuit,ref_m,100) setup=(circuit = copy(x_diag_circuit(1000)))


SUITE["pauliframe"]["quantumclifford"] = BenchmarkGroup(["quantumclifford"])
SUITE["pauliframe"]["quantumclifford"]["100"] = @benchmarkable  mctrajectory!(state, circuit) setup=(state= Register(one(Stabilizer, 101), [false]); circuit = copy(x_diag_circuit(100)))
SUITE["pauliframe"]["quantumclifford"]["1000"] = @benchmarkable mctrajectory!(state, circuit)  setup=(state= Register(one(Stabilizer, 1001), [false]); circuit = copy(x_diag_circuit(1000)))
