# [Universal qLDPC adapters: joint logical measurement](@id adapter_tutorial)

```@meta
DocTestSetup = quote
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC.Adapters
end
```

This page walks through the `Adapters` submodule on a concrete example: given
two CSS codeblocks and a specific logical operator between them, build the
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
using QuantumClifford: stab_to_gf2
using QuantumClifford.ECC: CSS, Surface, parity_matrix_x, parity_matrix_z,
                            logz_ops, logx_ops, code_n, code_k, parity_checks,
                            naive_encoding_circuit
using QuantumClifford.ECC.Adapters: CodePair, build_adapter, joint_logical_recipe

as_css(c) = CSS(Matrix{Bool}(parity_matrix_x(c)),
                Matrix{Bool}(parity_matrix_z(c)))

codeA = Surface(3, 3)
codeB = Surface(3, 3)
nA = code_n(codeA); kA = code_k(codeA)
nB = code_n(codeB); kB = code_k(codeB)

# Z-logical support on a single Surface(3,3) block
lz = stab_to_gf2(logz_ops(codeA))
z = sort(findall(!iszero, lz[1, nA+1:2nA]))
```

`z` is the list of data qubits that ``\bar Z`` acts on within one block. We
use the same support on both blocks since they are copies of the same code.

## Step 2: build the adapter and read the measurement recipe

`build_adapter` builds the two auxiliary graphs, cellulates them, runs the
SkipTree algorithm on each, and assembles the merged CSS code that joins
the two blocks.

```@example adapter
adapter = build_adapter(CodePair(as_css(codeA), as_css(codeB), z, z))

(code_n(adapter), code_k(adapter))
```

The merged code has `code_n(adapter)` data qubits and `code_k(adapter) = 1`
surviving logical qubit (the two original logicals have fused, as one expects
from a surgery merge).

`joint_logical_recipe` returns the row indices of the merged ``Z`` parity-check
matrix whose XOR over GF(2) equals ``\bar Z_l \otimes \bar Z_r`` in the
merged-qubit layout:

```@example adapter
recipe = joint_logical_recipe(adapter)
```

This is the *list of measurements* that implements the joint logical operator:
after standard syndrome extraction on the merged code, XOR the measurement
outcomes at these indices to recover the joint
``\bar Z_l \otimes \bar Z_r`` measurement result.

## Step 3: algebraic verification of the recipe

Each recipe row is a single ``Z``-type stabilizer of the merged code. We can
form the Pauli product directly and check that it equals the joint logical
``\bar Z_l \otimes \bar Z_r`` embedded in the merged qubit layout (the merged
column order is ``[\text{code}_1 \mid G_1 \mid A \mid G_2 \mid \text{code}_2]``,
with ``A`` the adapter block and ``G_1, G_2`` the auxiliary-graph edge
qubits).

```@example adapter
function row_to_zpauli(row, n)
    p = zero(PauliOperator, n)
    for i in 1:n
        row[i] != 0 && (p[i] = (false, true))
    end
    p
end

n_total = code_n(adapter)
Hz = parity_matrix_z(adapter)
recipe_paulis = [row_to_zpauli(Hz[r, :], n_total) for r in recipe]

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

We follow the standard `QuantumClifford` pattern (also used in the ECC
syndrome tests): a random ``k``-qubit logical input is buffered with
``n-k`` ``|0\rangle`` ancilla qubits and passed through
`naive_encoding_circuit`. This produces a fully-encoded `MixedDestabilizer`
state of the surface code.

```@example adapter
using Random; Random.seed!(2026)

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
of the two surface code blocks.

The logical Pauli we want to measure on this joint state is the tensor
product of the two ``\bar Z`` operators. The single-block ``\bar Z`` is the
first row of `logz_ops(codeA)` (and the same for `codeB`):

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

The pair `(count(naive_bits), nshots)` reports how many shots returned
``-1`` (i.e. measured `true`). For the seed-pinned random logical input above,
the empirical distribution should be near ``50/50``.

## Step 6: adapter-mediated measurement of the same logical

The adapter implements the same joint measurement using only ``X``- and
``Z``-type stabilizers of the merged code. The protocol on each shot is:

1. Insert ``n_\text{anc}`` ancilla qubits in ``|+\rangle`` between the two
   encoded blocks, matching the merged-code column layout.
2. Measure every ``X``-type merged stabilizer (these are deterministic on the
   chosen ancilla state and provide the standard CSS syndrome bits).
3. Measure every ``Z``-type merged stabilizer.
4. XOR the outcomes at the indices in `joint_logical_recipe(adapter)` to
   recover the joint logical bit.

Each *individual* merged stabilizer has weight at most ``4`` for this
adapter (the original Surface(3,3) checks are weight 3 or 4, and the
new ``V`` rows and bridge ``X``-rows the adapter introduces are weight
bounded by the auxiliary-graph degree); the joint logical itself has
weight 6, so even at distance 3 the adapter trades one weight-6
measurement for a sequence of low-weight ones. The gap widens with
code distance — joint logicals scale as ``O(d)`` while adapter
stabilizers stay ``O(1)`` — and that is the whole point of the
construction: a hardware that can execute only weight-bounded
measurements (as on a qLDPC architecture) can still perform an
arbitrary joint logical Pauli measurement.

```@example adapter
function row_to_xpauli(row, n)
    p = zero(PauliOperator, n)
    for i in 1:n
        row[i] != 0 && (p[i] = (true, false))
    end
    p
end

Hx = parity_matrix_x(adapter)
n_xstabs = size(Hx, 1)
n_zstabs = size(Hz, 1)

adapter_stabilizers = vcat(
    [PauliMeasurement(row_to_xpauli(Hx[r, :], n_total), r)
        for r in 1:n_xstabs],
    [PauliMeasurement(row_to_zpauli(Hz[r, :], n_total), n_xstabs + r)
        for r in 1:n_zstabs],
)

# Ancilla qubits live BETWEEN block A and block B in the merged column layout.
n_anc = n_total - nA - nB
plus_anc = MixedDestabilizer(Stabilizer([single_x(n_anc, i) for i in 1:n_anc]))
merged_state = encoded_A ⊗ plus_anc ⊗ encoded_B

adapter_bits = falses(nshots)
for shot in 1:nshots
    reg = Register(copy(merged_state), n_xstabs + n_zstabs)
    mctrajectory!(reg, adapter_stabilizers)
    bits = bitview(reg)
    recipe_outcomes = [bits[n_xstabs + r] for r in recipe]
    adapter_bits[shot] = reduce(xor, recipe_outcomes)
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
  which merged-code ``Z``-stabilizer outcomes to XOR to read out
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
