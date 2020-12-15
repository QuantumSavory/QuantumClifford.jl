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
# TODO
```

