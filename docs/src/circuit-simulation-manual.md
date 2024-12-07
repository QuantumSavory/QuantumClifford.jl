# [Circuit Simulation Manual](@id circuit-simulation-manual)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using StableRNGs
end
```

The library consists of two main parts: Tools for working with the algebra of Stabilizer tableaux and tools specifically for efficient Circuit Simulation. This chapter discusses the latter "higher level" circuit simulation tools.

# Classical Register

A [`Register`](@ref) encapsulates the state of a computer and includes both a tableaux and an array of classical bits, which can be used, for example, to store measurement outcomes. A [`MixedDestabilizer`](@ref) can be stored within the [`Register`](@ref), along with a set of classical bits for recording measurement results. The array of classical bits can be accessed through [`bitview`](@ref), while the tableaux can be examined using [`quantumstate`](@ref).

```jldoctest register
julia> s = MixedDestabilizer(T"YZ -XX XI IZ", 2)
𝒟ℯ𝓈𝓉𝒶𝒷
+ YZ
- XX
𝒮𝓉𝒶𝒷
+ X_
+ _Z

julia> reg = Register(s, [0,0])
Register{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}(MixedDestablizer 2×2, Bool[0, 0])

julia> quantumstate(reg)
𝒟ℯ𝓈𝓉𝒶𝒷
+ YZ
- XX
𝒮𝓉𝒶𝒷
+ X_
+ _Z

julia> stabilizerview(reg)
+ X_
+ _Z

julia> destabilizerview(reg)
+ YZ
- XX

julia> bitview(reg)
2-element Vector{Bool}:
 0
 0
```

Measurement results can be obtained using symbolic measurement operations such as [`sMX`](@ref), [`sMY`](@ref), and [`sMZ`](@ref), which can be applied with [`apply!`](@ref). 

```jldoctest
julia> rng = StableRNG(42); # hide

julia> md = MixedDestabilizer(T"XY -ZZ -XX -YZ", 2)
𝒟ℯ𝓈𝓉𝒶𝒷
+ XY
- ZZ
𝒮𝓉𝒶𝒷
- XX
- YZ

julia> reg = Register(md, [0, 0]);

julia> stabmx = apply!(rng, copy(reg), sMX(1, 1));

julia> quantumstate(stabmx)
𝒟ℯ𝓈𝓉𝒶𝒷
+ XY
- YZ
𝒮𝓉𝒶𝒷
- XX
- X_
```

Projective measurements with automatic phase randomization, including [`projectXrand!`](@ref), [`projectYrand!`](@ref), [`projectZrand!`](@ref), and [`projectrand!`](@ref) are available for the [`Register`](@ref) object.

```jldoctest
julia> rng = StableRNG(42); # hide

julia> s = MixedDestabilizer(T"YZ -XX XI IZ", 2)
𝒟ℯ𝓈𝓉𝒶𝒷
+ YZ
- XX
𝒮𝓉𝒶𝒷
+ X_
+ _Z

julia> reg = Register(s, [0, 0]);

julia> px = projectXrand!(rng, copy(reg), 2);

julia> quantumstate(px[1])
𝒟ℯ𝓈𝓉𝒶𝒷
+ Y_
+ _Z
𝒮𝓉𝒶𝒷
+ X_
- _X
```

# Pauli Frame

[`PauliFrame`](@ref) is a wrapper for a "frame" tableau. Each row represents the Pauli operation differing
the frame from the reference, behaving uniquely under Clifford operations.

```jldoctest frame
julia> circuit = [sHadamard(2), 
                  sHadamard(5), 
                  sCNOT(1, 2), 
                  sCNOT(2, 5), 
                  sMZ(1), 
                  sMZ(2), 
                  sMZ(4), 
                  sMZ(5)];

julia> frame = 5; qubits = 5; measurements=4;

julia> pframe = PauliFrame(frame, qubits, measurements);

julia> pframe = pftrajectories(pframe, circuit);

julia> pfmeasurements(pframe)
5×4 Matrix{Bool}:
 0  0  0  0
 0  0  0  0
 0  0  0  0
 0  0  0  0
 0  0  0  0

julia> sum(pframe.measurements)
0
```

Perturbative expansion with compactifying the same circuit used above with `compactify_circuit` and
`QuantumClifford._create_pauliframe` to create the [`PauliFrame`](@ref).

```jldoctest frame
julia> pfcircuit = compactify_circuit(circuit);

julia> pframe = QuantumClifford._create_pauliframe(pfcircuit; trajectories=5);

julia> result = pftrajectories(pframe, pfcircuit);

julia> pfmeasurements(result)
5×1 Matrix{Bool}:
 0
 0
 0
 0
 0

julia> sum(result.measurements)
0
```

Inject random Z errors on all frames and qubits with 0.5 probability.

```jldoctest 
julia> rng = StableRNG(42); # hide

julia> pframe = PauliFrame(4, 4, 2);

julia> QuantumClifford.initZ!(rng, pframe);

julia> pframe.frame
+ Z___
+ ___Z
+ __ZZ
+ ____
```

[`PauliFrame`](@ref) supports circuits with [`PauliError`](@ref), [`UnbiasedUncorrelatedNoise`](@ref), [`sMRZ`](@ref), [`sMRX`](@ref), [`sMX`](@ref), and [`ClassicalXOR`](@ref) for perturbative expansions. See [`applynoise!`](@ref), and [`pftrajectories`](@ref) as well.

For more examples, see the [notebook comparing the Monte Carlo and Perturbative method](https://nbviewer.jupyter.org/github/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Perturbative_Expansions_vs_Monte_Carlo_Simulations.ipynb) or this tutorial on [entanglement purification for many examples](https://github.com/QuantumSavory/QuantumClifford.jl/blob/master/docs/src/notebooks/Noisy_Circuits_Tutorial_with_Purification_Circuits.ipynb).
