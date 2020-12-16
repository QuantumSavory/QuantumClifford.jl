# Manual

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.Experimental.NoisyCircuits
end
```

Import with `using QuantumClifford.Experimental.NoisyCircuits`.

Here is an example of a purification circuit:

```@example
using QuantumClifford # hide
using QuantumClifford.Experimental.NoisyCircuits # hide
g1 = SparseGate(CNOT, [1,3]);
g2 = SparseGate(CNOT, [2,4]);
m = Measurement([X,X],[3,4]);
good_bell_state = S"XX
                           ZZ";
canonicalize_rref!(good_bell_state);
v = VerifyOp(good_bell_state,[1,2]);
n = NoiseOpAll(UnbiasedUncorrelatedNoise(0.01));
petrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v])
```

## Interface

`applyop_branches!(s::Stabilizer, g::Operation; max_order=1)::Vector{Tuple{Stabilizer,Int,Real,Int}}`
where the first `Int` is the status of the operation, the `Real` is the probability for that branch, and the second `Int` is the order of that branch.

There is also `applynoise_branches!` which is convenient for the `NoisyGate` and `NoisyMeasurement` and so on, but you can also just make up your own noise operator simply by implementing `applyop_branches!` for it.
