# Manual

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.Experimental.NoisyCircuits
end
```

Import with `using QuantumClifford.Experimental.NoisyCircuits`.

This module enables the simulation of noisy Clifford circuits through a perturbative expansion in the noise parameter (assuming the noise is small). Instead of simulating many Monte Carlo trajectories, only the leading order trajectories are exhaustively enumerated and simulated.

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

petrajectories(initial_state, circuit)
```

## Interface

`applyop_branches!(s::Stabilizer, g::Operation; max_order=1)::Vector{Tuple{Stabilizer,Int,Real,Int}}`
where the first `Int` is the status of the operation, the `Real` is the probability for that branch, and the second `Int` is the order of that branch.

There is also `applynoise_branches!` which is convenient for the `NoisyGate` and `NoisyMeasurement` and so on, but you can also just make up your own noise operator simply by implementing `applyop_branches!` for it.
