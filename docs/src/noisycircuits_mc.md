# [Monte Carlo simulations of noisy Clifford circuits](@id noisycircuits_mc)

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

This module enables the simulation of noisy Clifford circuits through a Monte Carlo method where the same circuit is evaluated multiple times with random errors interspersed through it as prescribed by a given error model.

Below is an example of a purification circuit. We first prepare the circuit we desire to use, including a noise model. `Quantikz.jl` was is used to visualize the circuit.

```@example 1
using QuantumClifford # hide
using QuantumClifford.Experimental.NoisyCircuits # hide
using Quantikz # hide
good_bell_state = S"XX
                    ZZ"
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
```

And we can run a Monte Carlo simulation of that circuit with [`mctrajectories`](@ref).

```@example 1
mctrajectories(initial_state, circuit, trajectories=500)
```
... and without purification

```@example 1
mctrajectories(initial_state, [n, v], trajectories=500)
```

... and without noise

```@example 1
mctrajectories(initial_state, [g1,g2,m,v], trajectories=500)
```

... and use of [`Register`](@ref) and [`SparseGate`](@gate) for Monte Carlo simulation of that circuit with [`mctrajectories`](@ref).

```@example 2
g1 = SparseGate(tCNOT, [1,3])
g2 = SparseGate(tCNOT, [2,4])
m = BellMeasurement([sMX(3),sMX(4)])
good_bell_state = ghz(2)
canonicalize_rref!(good_bell_state)
v = VerifyOp(good_bell_state,[1,2])
n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.03))
init = Register(MixedDestabilizer(good_bell_state⊗good_bell_state))
with_purification = mctrajectories(init, [n,g1,g2,m,v], trajectories=500)
```

For more examples, see the [notebook comparing the Monte Carlo and Perturbative method](https://nbviewer.jupyter.org/github/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Perturbative_Expansions_vs_Monte_Carlo_Simulations.ipynb) or this tutorial on [entanglement purification for many examples](https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Noisy_Circuits_Tutorial_with_Purification_Circuits.ipynb).

## Interface for custom operations

If you want to create a custom gate type (e.g. calling it `Operation`), you need to definite the following methods.

`applywstatus!(s::T, g::Operation)::Tuple{T,Symbol}` where `T` is a tableaux type like [`Stabilizer`](@ref) or a [`Register`](@ref).
The `Symbol` is the status of the operation. Predefined statuses are kept in the `registered_statuses` list, but you can add more.
Be sure to expand this list if you want the trajectory simulators using your custom statuses to output all trajectories.

There is also [`applynoise!`](@ref) which is convenient wait to create a noise model that can then be plugged into the [`NoisyGate`](@ref) struct,
letting you reuse the predefined perfect gates and measurements.
However, you can also just make up your own noise operator simply by implementing [`applywstatus!`](@ref) for it.

You can also consult the [list of implemented operators](@ref noisycircuits_ops).