# Simulation of Noisy Clifford Circuits

```@meta
DocTestSetup = quote
    using QuantumClifford
end
```

We have support for the simulation of noisy Clifford circuits, providing a framework for exploring noise processes within quantum systems

Both [Monte Carlo](@ref noisycircuits_mc) and [Perturbative Expansion](@ref noisycircuits_perturb) approaches are supported. When performing a perturbative expansion in the noise parameter, the expansion can optionally be performed symbolically, to arbitrary high orders.

Multiple [notebooks with examples](@ref tutandpub) are also available.
For instance, see this tutorial on [entanglement purification for many examples](https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Noisy_Circuits_Tutorial_with_Purification_Circuits.ipynb).
