# [Noise Processes](@id noise)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using StableRNGs
    rng = StableRNG(42)
end
```

As seen in the list of [possible gates](@ref all-operations),
the simulator is capable of modeling different types of noise.
If that is your goal, please consider using the available
[Monte Carlo simulator](@ref noisycircuits_mc) or the
[Symbolic Perturbative Expansion system](@ref noisycircuits_perturb).

The implemented types of noise include:

- [`UnbiasedUncorrelatedNoise`](@ref)
- [`PauliNoise`](@ref)

The low-level functionality to work with noise is `applynoise!`,
but most of the time you would probably just want to use
[`PauliError`](@ref),  [`NoisyGate`](@ref), [`NoiseOp`](@ref) and [`NoiseOpAll`](@ref).