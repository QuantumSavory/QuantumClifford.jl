# Simulation of Noisy Clifford Circuits

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.Experimental.NoisyCircuits
end
```

!!! warning "Unstable"
    This is unfinished experimental functionality that will change significantly.

We have experimental support for simulation of noisy Clifford circuits which can be imported with `using QuantumClifford.Experimental.NoisyCircuits`.

Both [Monte Carlo](@ref noisycircuits_mc) and [Perturbative Expansion](@ref noisycircuits_perturb) approaches are supported. When performing a perturbative expansion in the noise parameter, the expansion can optionally be performed symbolically, to arbitrary high orders.

Multiple [notebooks with examples](@ref tutandpub) are also available.
For instance, see this tutorial on [entanglement purification for many examples](https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Noisy_Circuits_Tutorial_with_Purification_Circuits.ipynb).