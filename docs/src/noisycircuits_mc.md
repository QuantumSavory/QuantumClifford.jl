# Monte Carlo simulations of noisy Clifford circuits

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.Experimental.NoisyCircuits
end
```

Import with `using QuantumClifford.Experimental.NoisyCircuits`.

This module enables the simulation of noisy Clifford circuits through a Monte Carlo method where the same circuit is evaluated multiple times with random errors interspersed through it as prescribed by a given error model.

Here is an example of a purification circuit:

```@example
using QuantumClifford # hide
using QuantumClifford.Experimental.NoisyCircuits # hide
good_bell_state = S"XX
                    ZZ"
canonicalize_rref!(good_bell_state)
initial_state = good_bell_state⊗good_bell_state

g1 = SparseGate(CNOT, [1,3]) # CNOT between qubit 1 and qubit 3 (both with Alice)
g2 = SparseGate(CNOT, [2,4]) # CNOT between qubit 2 and qubit 4 (both with Bob)
m = BellMeasurement([X,X],[3,4]) # Bell measurement on qubit 3 and 4
v = VerifyOp(good_bell_state,[1,2]) # Verify that qubit 1 and 2 indeed form a good Bell pair
epsilon = 0.01 # The error rate
n = NoiseOpAll(UnbiasedUncorrelatedNoise(epsilon))

# This circuit performs a depolarization at rate `epsilon` to all qubits,
# then bilater CNOT operations
# then a Bell measurement
# followed by checking whether the final result indeed corresponds to the correct Bell pair.
circuit = [n,g1,g2,m,v]

mctrajectories(initial_state, circuit, trajectories=500)
```

For more examples, see the [notebook comparing the Monte Carlo and Perturbative method](https://nbviewer.jupyter.org/github/Krastanov/QuantumClifford.jl/blob/master/docs/src/notebooks/Perturbative_Expansions_vs_Monte_Carlo_Simulations.ipynb).

## Interface

`applyop!(s::Stabilizer, g::Operation)::Tuple{Stabilizer,Int}`
where the `Int` is the status of the operation. Predefined statuses are kept in the `statuses` dictionary:
```julia
const statuses = Dict(0=>:continue,
                      1=>:detected_failure,
                      2=>:undetected_failure,
                      3=>:true_success)
const s_continue = 0
const s_detected_failure = 1
const s_undetected_failure = 2
const s_true_success = 3
```

Be sure to expand this dictionary if you want the trajectory simulators using your custom statuses to output all trajectories.

TODO: We probably do not want to use this dictionary, rather a type/dispatch-based approach or symbols/namedtuple-based approach? Something to consider for a future version.

There is also `applynoise!` which is convenient wait to create a noise model that can then be plugged into the `NoisyGate`, `NoisyMeasurement`, etc, letting you reuse the predefined perfect gates and measurements. However, you can also just make up your own noise operator simply by implementing `applyop!` for it.