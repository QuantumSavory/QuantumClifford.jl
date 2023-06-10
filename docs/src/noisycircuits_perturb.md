# [Perturbative expansions for simulating noisy Clifford circuits](@id noisycircuits_perturb)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.Experimental.NoisyCircuits
    using Quantikz
end
CurrentModule = QuantumClifford.Experimental.NoisyCircuits
```

!!! warning "Unstable"
    This is experimental functionality with an unstable API.
    
Import with `using QuantumClifford.Experimental.NoisyCircuits`.

This module enables the simulation of noisy Clifford circuits through a perturbative expansion in the noise parameter (assuming the noise is small).
Instead of simulating many Monte Carlo trajectories, only the leading order trajectories are exhaustively enumerated and simulated.

Here is an example of a purification circuit (the same circuit seen in the [Monte Carlo example](@ref noisycircuits_mc))

```@example
using QuantumClifford # hide
using QuantumClifford.Experimental.NoisyCircuits # hide
using Quantikz # hide
good_bell_state = S"XX
                    ZZ"
canonicalize_rref!(good_bell_state)
initial_state = MixedDestabilizer(good_bell_state⊗good_bell_state)

g1 = sCNOT(1,3) # CNOT between qubit 1 and qubit 3 (both with Alice)
g2 = sCNOT(2,4) # CNOT between qubit 2 and qubit 4 (both with Bob)
m = BellMeasurement([sMX(3),sMX(4)]) # Bell measurement on qubit 3 and 4
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

For more examples, see the [notebook comparing the Monte Carlo and Perturbative method](https://nbviewer.jupyter.org/github/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Perturbative_Expansions_vs_Monte_Carlo_Simulations.ipynb) or this tutorial on [entanglement purification](https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Noisy_Circuits_Tutorial_with_Purification_Circuits.ipynb).

## Symbolic expansions

The perturbative expansion method works with symbolic variables as well. One can use any of the symbolic libraries available in Julia and simply plug symbolic parameters in lieu of numeric parameters. A detailed example is available as a [Jupyter notebook](https://nbviewer.jupyter.org/github/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Symbolic_Perturbative_Expansions.ipynb).

## Interface for custom operations

If you want to create a custom gate type (e.g. calling it `Operation`), you need to definite the following methods.

`applyop_branches!(s::T, g::Operation; max_order=1)::Vector{Tuple{T,Symbol,Real,Int}}` where `T` is a tableaux type like [`Stabilizer`](@ref) or a [`Register`](@ref).
The `Symbol` is the status of the operation, the `Real` is the probability for that branch, and the `Int` is the order of that branch.

There is also `applynoise_branches!` which is convenient for use in `NoisyGate`, but you can also just make up your own noise operator simply by implementing `applyop_branches!` for it.

You can also consult the [list of implemented operators](@ref noisycircuits_ops).