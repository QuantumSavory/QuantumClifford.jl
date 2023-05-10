# [Pauli Frames](@id pauliframe)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using Random
end
```

A Pauli Frame simulation can be created by interacting directly with the struct and its methods, or by calling the 
[`pauliFrameCircuitHandler`](@ref) or [`mctrajectory!`](@ref) function, which will run an entire circuit for you, calling the appropiate methods and in the case of [`pauliFrameCircuitHandler`](@ref), it will create and initialize a [`PauliFrame`](@ref) object. If you do create your own [`PauliFrame`](@ref) object, please consider calling [`initZ!`](@ref) on it so that it can properly sample non-deterministic circuits.

Currently, all 1 and 2 qubit gates from QuantumClifford work. However, only Z basis measurement is available.
Another current assumption about measurement is that all measurements happen at the end of the circuit. 

The purpose of a Pauli frame simulation is to reduce sample times for a circuit, at the cost of having to provide a set
of reference measurements that were obtained by another means. Therefore, for this to be of any use, it has to be faster
than the alternative, which in the case of QuantumClifford.jl, would be [`mctrajectory`](@ref). Using `pauli_frame_benchmarks.jl`, we can see that while simulating an entire 1,000 qubit 2,000 gate circuit without Pauli frames takes about 12.042 ms and 500.78 KiB, using Pauli frames, one can generate 1,000 samples of the very same circuit with only 5.402 ms and 726.00 KiB!