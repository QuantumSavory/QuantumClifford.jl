# Cold-start benchmark

This harness compares package-cache build, import, and first-use latency. It is
separate from the BenchmarkTools and AirspeedVelocity suite, which measures
steady-state operations after QuantumClifford is loaded.

From the repository root, give the harness a new output directory plus two or
more `label=checkout` variants:

```sh
benchmark/precompile/run.sh /tmp/quantumclifford-precompile-results \
    base=/path/to/base/QuantumClifford.jl \
    head=/path/to/head/QuantumClifford.jl
```

The first variant is the default baseline. Each variant must be a clean,
committed checkout, and all variants must have identical `Project.toml` files
apart from the top-level package version.
The output directory must be outside every measured checkout. The harness
refuses to overwrite existing result files. By default, it resolves one
consumer Manifest and points it at each checkout in turn. It creates one seed
depot with dependency caches and a new writable depot for every QuantumClifford
package-cache build. Dependency setup may access package servers. After setup,
cache builds and samples run with package offline mode enabled. For reportable
runs, seed setup loads QuantumClifford from a detached Git archive with separate
source-file inodes, then deletes that package cache. It therefore does not
page-warm a measured checkout. A non-reportable run with a dirty baseline uses
a detached copy of that fixed working state so setup and measurement still use
the same source.
Before timing, every measured variant gets one discarded QuantumClifford cache
build in a fresh overlay depot through the same stable checkout link. Each
discarded build must emit cache bytes. The harness orders these builds by a
SHA-256 value derived from the commit and initial checkout-state hash,
independent of baseline/candidate role. Identical source states retain argument
order because their filesystem exposure is equivalent.

Variant labels and scenario names must start with an ASCII letter or digit and
contain only ASCII letters, digits, dots, underscores, or hyphens. This keeps
the comma-separated controls, TSV output, metadata keys, and Markdown output
unambiguous.

Use environment variables to select the amount of work and the comma-separated
scenario list:

```sh
QC_PRECOMPILE_BUILDS=5 \
QC_PRECOMPILE_SAMPLES=4 \
QC_PRECOMPILE_SCENARIOS=tableau,ecc \
benchmark/precompile/run.sh /tmp/quantumclifford-precompile-results \
    base=/path/to/base head=/path/to/head
```

Use five cache builds and four recorded samples for reportable candidate
measurements. The default of one build and two samples is intended only for
smoke tests. Supported scenarios are `pauli`, `tableau`, and `ecc`.

For a multi-candidate experiment, use `QC_PRECOMPILE_EXTRA_SCENARIOS` to run a
motivating scenario only for its named variant and for the baseline. This keeps
the common headline scenarios on every variant:

```sh
QC_PRECOMPILE_SCENARIOS=tableau,ecc \
QC_PRECOMPILE_EXTRA_SCENARIOS=pauli_only=pauli \
benchmark/precompile/run.sh /tmp/quantumclifford-precompile-results \
    base=/path/to/base pauli_only=/path/to/pauli-candidate
```

By default, the first variant is the baseline for every candidate. Use
`QC_PRECOMPILE_BASELINES` when candidates need different baselines, such as a
cumulative experiment in which every stage is compared with the preceding
stage:

```sh
QC_PRECOMPILE_BASELINES=stage2=stage1,stage3=stage2 \
benchmark/precompile/run.sh /tmp/quantumclifford-precompile-results \
    base=/path/to/base \
    stage1=/path/to/stage1 \
    stage2=/path/to/stage2 \
    stage3=/path/to/stage3
```

To keep exact dependency versions across separate harness invocations, reuse
the normalized consumer files from the first result directory:

```sh
QC_PRECOMPILE_CONSUMER_PROJECT=/path/to/first-results/consumer-Project.toml \
QC_PRECOMPILE_CONSUMER_MANIFEST=/path/to/first-results/consumer-Manifest.toml \
benchmark/precompile/run.sh /tmp/quantumclifford-precompile-later-results \
    base=/path/to/base \
    head=/path/to/head
```

Set both variables together and use a new output directory. The inputs must be
distinct regular files. The Project must match the harness-generated
QuantumClifford consumer Project. The normalized
Manifest must contain exactly one `__QUANTUMCLIFFORD_CHECKOUT__` path in its
QuantumClifford entry. Setup materializes that placeholder as the stable
temporary checkout link. It first uses `Pkg.is_manifest_current` in offline
mode to require a Manifest resolved from that Project, then runs only
`Pkg.instantiate()` against the saved dependency graph. It does not resolve or
update dependencies. The copied consumer files in the new result directory
must remain byte-identical to the inputs or the harness fails.

The default scenarios are `tableau` and `ecc`. The first covers core tableau
algebra and projection; the second covers Steane-code circuits, Pauli-frame
simulation, and Shor-code table decoding. Each scenario also gets one
discarded filesystem warm-up per build. The harness fixes Julia,
package-precompile, BLAS, and OpenMP thread counts to one; disables startup and
history files; uses `JULIA_LOAD_PATH=@:@stdlib`; and clears inherited Julia CPU
target, project, and depot overrides. It requires GNU/Linux and Julia 1.12.6.
Set `JULIA` to select that Julia executable. A different Julia
version or a dirty checkout is allowed only for a non-reportable smoke run by
setting `QC_PRECOMPILE_ALLOW_JULIA_MISMATCH=1` or
`QC_PRECOMPILE_ALLOW_DIRTY=1`, respectively. Enabling either override records
`reportable=false` and a reason in `metadata.txt`, even when the current Julia
version matches or the checkout is clean. A run also records
`reportable=false` when it uses fewer than five builds, fewer than four samples
per build, or omits either `tableau` or `ecc` from the common scenarios.
Thus the small default run and the descriptive pull-request workflow are
explicitly non-reportable. A run that satisfies these repetition and headline
requirements without an override records `reportable=true` and an empty
`nonreportable_reasons` list.

The harness commit and its three source hashes are checked before setup and
again before summarization. Recorded processes execute temporary read-only
snapshots of `scenarios.jl` and `summarize.jl`; their hashes are
checked after creation and before summarization. An edit in the harness
checkout therefore cannot mix source versions within a long run. The harness
also snapshots a Git worktree-state hash for every measured checkout, including
tracked differences and nonignored untracked contents. It rechecks the hash
before each cache build and after all measurements. This rejects source mutation
during an allow-dirty run; the override permits a fixed initial dirty state only.

Each candidate gets a separate comparison, with nearby independent cache
builds for the candidate and its mapped baseline. To counterbalance systematic
drift, odd-numbered builds visit candidate arguments in ascending order and
even-numbered builds visit them in descending order. Candidate indices are
their one-based positions after the baseline argument. Within a comparison,
the baseline runs before the candidate when `build + candidate_index` is odd;
the candidate runs first when it is even. Five-build reportable runs therefore
use a deterministic 3/2 split of pair order, while preserving adjacency. The
metadata records these policies, each candidate index, and the exact comparison
and pair sequence for every build as `schedule.build.N`.

The output directory contains `raw.tsv`, the per-build `build-summary.tsv`, the
aggregate `summary.tsv`, the rendered `summary.md`, `metadata.txt`, the saved
consumer Project `consumer-Project.toml`, and the normalized consumer Manifest
`consumer-Manifest.toml`. Each `build-summary.tsv` row includes the raw
`build_seconds` and integer `cache_bytes` values for that build, in addition to
the scenario medians. The aggregate
summary reports medians, interquartile ranges, and the number of builds with a
material total-latency change. A material change is at least 50 ms and 5% is
used when it is larger. Performance differences are descriptive; scenario
assertion, package-cache, dependency-control, or harness failures return a
nonzero status.

The total-latency timer starts immediately before `using QuantumClifford` and
stops after the first scenario call. It therefore includes the small amount of
harness setup between the separately reported import and first-call timers.
The `total_metric` metadata field records this definition.

The copied Manifest uses `__QUANTUMCLIFFORD_CHECKOUT__` as a path placeholder;
replace it with the checkout used for reproduction. The `manifest_sha256`
metadata field hashes this normalized copy. The metadata also records the
consumer-environment mode, consumer Project hash, pre-normalization Manifest
hash, harness commit, harness file hashes, reportability state, checkout state
hashes, and the discarded cache-warm-up count, policy, and order. Reuse mode
additionally records the canonical source paths and hashes for both inputs. The pre-normalization
`resolved_manifest_sha256` includes an invocation-specific temporary checkout
path and is expected to change across reuse invocations; use `manifest_sha256`
as the stable normalized dependency identity.

See the [Julia command-line reference](https://docs.julialang.org/en/v1/manual/command-line-interface/)
for the compilation controls and the [PrecompileTools workload guide](https://julialang.github.io/PrecompileTools.jl/stable/)
for workload design.
