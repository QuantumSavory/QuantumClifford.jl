# Monte Carlo simulations of noisy circuits

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
mctrajectories(good_bell_state⊗good_bell_state, [n,g1,g2,m,v], trajectories=500)
```

## Interface

`applyop!(s::Stabilizer, g::Operation)::Tuple{Stabilizer,Int}`
where the `Int` is the status of the operation. Predefined statuses are kept in the `statuses` dictionary:
```julia
const statuses = Dict(0=>:continue,
                      1=>:detected_failure,
                      2=>:undetected_failure,
                      3=>:true_success)
const s_continue = 0
const s_detected_failure = 1
const s_undetected_failure = 2
const s_true_success = 3
```

Be sure to expand this dictionary if you want the trajectory simulators using your custom statuses to output all trajectories.

TODO: We probably do not want to use this dictionary, rather a type/dispatch-based approach or symbols/namedtuple-based approach? Something to consider for a future version.

There is also `applynoise!` which is convenient wait to create a noise model that can then be plugged into the `NoisyGate`, `NoisyMeasurement`, etc, letting you reuse the predefined perfect gates and measurements. However, you can also just make up your own noise operator simply by implementing `applyop!` for it.