# [Universal qLDPC adapters: joint logical measurement](@id adapter_tutorial)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC.Adapters
end
```

This page walks through the `Adapters` submodule on a concrete example: given
two CSS codeblocks and a specific logical operator on each, build the
universal adapter of [swaroop2026universal](@cite), read off the list of
stabilizer measurements that implements the joint logical operator, and
compare against the naive direct-measurement protocol in simulation.

The construction realises a **joint Pauli measurement** ``\bar Z_l \otimes \bar Z_r``
between a logical of one codeblock and a logical of another, using only
weight-bounded ``X``- and ``Z``-type stabilizer measurements on a merged code.
This is one half of a code-surgery protocol — the other half (the *split*
step that restores two independent codeblocks afterwards) is not yet
implemented and is called out at the end.

## Step 1: pick the two codes and the logicals to fuse

We use two independent copies of the ``3 \times 3`` planar surface code
(`Surface(3, 3)` from `QECCore`) — `code_n = 13`, `code_k = 1`,
`distance = 3` — and target the ``\bar Z`` logical on each block.

```@example adapter
using QuantumClifford
using QuantumClifford.ECC

codeA = Surface(3, 3)
codeB = Surface(3, 3)
nA = code_n(codeA); kA = code_k(codeA)
nB = code_n(codeB); kB = code_k(codeB)

# Z-logical operator to adapt
z = logz_ops(codeA)[1]
```

## Step 2: build the adapter and read the measurement recipe

`build_adapter` builds the two auxiliary graphs, cellulates them, runs the
SkipTree algorithm on each, and assembles the merged CSS code that joins
the two blocks. The two codes and the two logical operators are passed
directly — no manual conversion required.

```@example adapter
adapter = build_adapter(codeA, codeB, z, z)

(code_n(adapter), code_k(adapter))
```

The merged code has `code_n(adapter)` data qubits and `code_k(adapter) = 1`
surviving logical qubit (the two original logicals have fused, as one expects
from a surgery merge).

`joint_logical_recipe` returns the row indices of the merged code's
parity-check tableau whose product is the joint logical
``\bar Z_l \otimes \bar Z_r``:

```@example adapter
recipe = joint_logical_recipe(adapter)
```

The indices correspond to rows in `parity_checks(adapter)`. After standard
syndrome extraction on the merged code, XOR the measurement outcomes at
these indices to recover the joint ``\bar Z_l \otimes \bar Z_r``
measurement result.

## Step 3: algebraic verification of the recipe

This step is simply a verification for correctness, not necessary in actual use.

This step is simply a verification for correctness — not necessary in actual
use. Each recipe row is a single ``Z``-type stabilizer of the merged code;
their Pauli product equals the joint logical
``\bar Z_l \otimes \bar Z_r`` embedded in the merged qubit layout (the
column order is ``[\text{code}_1 \mid G_1 \mid A \mid G_2 \mid \text{code}_2]``,
with ``A`` the adapter block and ``G_1, G_2`` the auxiliary-graph edge
qubits).

```@example adapter
H = parity_checks(adapter)
recipe_paulis = H[recipe]

prod(recipe_paulis)
```

The product is supported only on data qubits 3, 6, 9 (the ``\bar Z`` support
of block 1) and 23, 26, 29 (the ``\bar Z`` support of block 2 — offset by
``n_A + n_\text{anc} = 13 + 7 = 20``). All auxiliary qubits are trivial. This
is exactly ``\bar Z_l \otimes \bar Z_r``.

The cancellation on the auxiliary qubits is by design: the ``V_l`` row block
contributes ``Z`` on every column of the adapter sub-block ``A`` (because the
SkipTree permutation ``P_l`` places exactly one ``1`` per ``A`` column), and
the ``V_r`` row block contributes the same ``Z`` on ``A``. Since
``Z^2 = I`` over ``\mathbb F_2``, the auxiliary support cancels and only the
joint logical survives on the data qubits.

## Step 4: prepare encoded logical states in simulation

The next two steps are also a verification of correctness — this time by
simulating either the measurement of the adapter we have built, or the much
more naive (and not fault tolerant or easy to implement on real hardware)
direct measurement of the large cross-code logical operator without the use
of an adapter.

We follow the standard `QuantumClifford` pattern (also used in the ECC
syndrome tests): a random ``k``-qubit logical input is buffered with
``n-k`` ``|0\rangle`` ancilla qubits and passed through
`naive_encoding_circuit`. This produces a fully-encoded `MixedDestabilizer`
state of the surface code.

```@example adapter
using Random; Random.seed!(2026) # hide

physical_state_A = one(Stabilizer, nA - kA) ⊗ random_stabilizer(kA)
physical_state_B = one(Stabilizer, nB - kB) ⊗ random_stabilizer(kB)

encoded_A = MixedDestabilizer(physical_state_A)
mctrajectory!(encoded_A, naive_encoding_circuit(codeA))

encoded_B = MixedDestabilizer(physical_state_B)
mctrajectory!(encoded_B, naive_encoding_circuit(codeB))

full_state = encoded_A ⊗ encoded_B
nqubits(full_state)
```

`full_state` is the joint encoded state on the ``n_A + n_B = 26`` data qubits
of the two surface code blocks. The logical Pauli we want to measure on this
joint state is the tensor product of the two ``\bar Z`` operators:

```@example adapter
logicalA = logz_ops(codeA)[1]
logicalB = logz_ops(codeB)[1]
full_logical = logicalA ⊗ logicalB
```

## Step 5: naive measurement of the joint logical

The *naive* protocol is to construct ``\bar Z_l \otimes \bar Z_r`` as a single
weight-6 Pauli and measure it directly with [`PauliMeasurement`](@ref) on a
[`Register`](@ref) wrapping `full_state`. We harvest one classical bit per
trajectory.

```@example adapter
nshots = 1000

naive_bits = falses(nshots)
for shot in 1:nshots
    reg = Register(copy(full_state), 1)
    mctrajectory!(reg, [PauliMeasurement(full_logical, 1)])
    naive_bits[shot] = bitview(reg)[1]
end

(count(naive_bits), nshots)
```

The pair `(count(naive_bits), nshots)` reports how many shots measured
`true`. The fraction will be 0 or 1 for a deterministic measurement, or
near 0.5 for a logical operator that does not stabilize the current
logical state.

## Step 6: adapter-mediated measurement of the same logical

Here we get back to actually using the adapter, the proper way to do
fault-tolerant logical measurements, and see that it gives the same result
as the naive method from above. This adapter approach is what you would
actually want when doing real QEC modelling work.

The adapter implements the same joint measurement using only ``X``- and
``Z``-type stabilizers of the merged code. The protocol on each shot is:

1. Insert ``n_\text{anc}`` ancilla qubits in ``|+\rangle`` between the two
   encoded blocks, matching the merged-code column layout.
2. Measure every merged stabilizer (the ``X``-stabilizers are deterministic
   on the chosen ancilla state and provide the standard CSS syndrome bits;
   the ``Z``-stabilizers include the joint-logical recipe).
3. XOR the outcomes at the indices in `joint_logical_recipe(adapter)` to
   recover the joint logical bit.

Each *individual* merged stabilizer has weight at most ``4`` for this
adapter; the joint logical itself has weight ``6``, so even at distance
``3`` the adapter trades one weight-``6`` measurement for a sequence of
low-weight ones. The gap widens with code distance — joint logicals scale
as ``O(d)`` while adapter stabilizers stay ``O(1)`` — and that is the whole
point of the construction: hardware that can execute only weight-bounded
measurements (as on a qLDPC architecture) can still perform an arbitrary
joint logical Pauli measurement.

```@example adapter
n_total = code_n(adapter)
n_anc = n_total - nA - nB
plus_anc = MixedDestabilizer(Stabilizer([single_x(n_anc, i) for i in 1:n_anc]))
merged_state = encoded_A ⊗ plus_anc ⊗ encoded_B

adapter_stabilizers = [PauliMeasurement(H[r], r) for r in eachindex(H)]

adapter_bits = falses(nshots)
for shot in 1:nshots
    reg = Register(copy(merged_state), length(H))
    mctrajectory!(reg, adapter_stabilizers)
    bits = bitview(reg)
    adapter_bits[shot] = reduce(xor, bits[r] for r in recipe)
end

(count(adapter_bits), nshots)
```

## Step 7: compare the two protocols

The joint-logical bit distributions sampled by the two protocols agree
within Monte Carlo error:

```@example adapter
(naive_p1 = count(naive_bits) / nshots,
 adapter_p1 = count(adapter_bits) / nshots)
```

Per-shot the outcomes do not agree (the two protocols sample independent
copies of `full_state`), but the empirical Bernoulli parameter is the same:
the two protocols implement the same logical measurement.

On an eigenstate of the joint logical the equivalence is sharper. If we
instead prepare the encoded logical ``|0\rangle \otimes |0\rangle`` — a
``+1`` eigenstate of ``\bar Z_l \otimes \bar Z_r`` — both protocols return
``0`` on every shot; if we then apply ``\bar X_l`` on block ``A``, both
return ``1`` on every shot. This deterministic agreement follows from the
algebraic recipe identity in Step 3: the XOR of the recipe rows is the
joint logical Pauli up to merged-code stabilizers, so on an eigenstate the
two protocols are forced to report the same bit.

## Caveat: the splitting step is not yet implemented

The full fault-tolerant code-surgery protocol summarised at the end of
Section II of [swaroop2026universal](@cite) (the auxiliary-graph surgery
protocol of Williamson–Yoder, Ref. [24] in that paper) has four steps:

1. **Initialise** all ancilla edge qubits in ``|+\rangle``.
2. **Repeat** the measurement of the merged code's stabilizers at least
   ``d`` times (for phenomenological fault distance ``d``).
3. **Split**: measure out all the edge qubits in the ``X`` basis.
4. **Correct**: apply a Pauli correction on the original code qubits to
   return to the original codespace.

What this PR provides:

- The merged CSS code itself: `build_adapter` returns an
  `AbstractCSSCode` whose `parity_matrix_x` / `parity_matrix_z` /
  `code_n` / `code_s` plug directly into the rest of the QuantumClifford
  ECC infrastructure, so existing decoders (e.g. the BP-OSD harness via
  `PyQDecoders`) operate on it without an adapter shim.
- The joint-logical recipe: `joint_logical_recipe(adapter)` identifies
  which merged-code stabilizer outcomes to XOR to read out
  ``\bar Z_l \otimes \bar Z_r``.
- The simulation example above, which is a single-round noiseless
  realisation of steps 1 and 2.

What is **not yet implemented**: the explicit ``d``-round repeated
measurement schedule, the ``X``-basis edge-qubit readout that performs
the split (step 3), and the Pauli frame correction (step 4). These are
planned for a follow-up PR. Without the split, the example above can
only be run as a one-shot joint measurement followed by a fresh state
preparation; combining adapter measurements with mid-circuit logical
operations across multiple rounds will require the split routine.
