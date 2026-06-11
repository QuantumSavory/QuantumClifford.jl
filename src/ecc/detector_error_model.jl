"""
$TYPEDEF

An in-memory representation of a code-capacity detector error model, a list of independent
single-qubit Pauli error mechanisms together with the detectors (parity checks) and logical
observables that each mechanism flips. It can be serialized to Stim's detector error model
(`.dem`) text format with [`write_detector_error_model`](@ref), e.g. for use with external
Stim-compatible decoders.

Create it from an ECC code with [`detector_error_model`](@ref).

### Fields

$TYPEDFIELDS

See also: [`detector_error_model`](@ref), [`write_detector_error_model`](@ref), [`faults_matrix`](@ref)
"""
struct DetectorErrorModel
    """The number of detectors, equal to the number of rows of `parity_checks(code)`. Detector `Dᵢ` (0-based, as in Stim) corresponds to row `i+1` (1-based) of the parity check tableau."""
    num_detectors::Int
    """The number of logical observables, equal to the number of rows of `faults_matrix(code)` (twice the number of logical qubits). Observable `Lⱼ` (0-based) corresponds to row `j+1` (1-based) of the faults matrix, i.e. `L₀ … L_{k-1}` are the logical X observables and `L_k … L_{2k-1}` are the logical Z observables of a code with `k` logical qubits."""
    num_observables::Int
    """The error mechanisms: one entry per enabled (nonzero-probability) single-qubit Pauli fault, each a named tuple of the fault `probability`, the 0-based indices of the flipped `detectors`, and the 0-based indices of the flipped `observables`. Mechanisms are ordered by qubit index, then by Pauli kind in `X`, `Y`, `Z` order; the index lists are sorted ascending. Mechanisms that flip no detector and no observable are omitted."""
    errors::Vector{@NamedTuple{probability::Float64, detectors::Vector{Int}, observables::Vector{Int}}}
end

function Base.:(==)(a::DetectorErrorModel, b::DetectorErrorModel)
    a.num_detectors == b.num_detectors && a.num_observables == b.num_observables && a.errors == b.errors
end
Base.hash(dem::DetectorErrorModel, h::UInt) = hash((dem.num_detectors, dem.num_observables, dem.errors), hash(DetectorErrorModel, h))

function Base.show(io::IO, ::MIME"text/plain", dem::DetectorErrorModel)
    print(io, "DetectorErrorModel with $(dem.num_detectors) detectors, $(dem.num_observables) logical observables, and $(length(dem.errors)) error mechanisms")
end

"""
$TYPEDSIGNATURES

Compute the code-capacity detector error model of a code, suitable for export to Stim's
`.dem` format via [`write_detector_error_model`](@ref).

The model contains one independent error mechanism per physical qubit per enabled Pauli
fault: a fault is enabled by giving it a nonzero probability through the keyword arguments
`px`, `py`, `pz` (each a probability in [0, 1], all zero by default). Each mechanism is an
independent Bernoulli event, matching the semantics of Stim's `error(p)` instruction --
in particular `px`, `py`, `pz` are probabilities of three *independent* faults on each
qubit, not the branch probabilities of a single-qubit depolarizing channel.

For a fault `P` on qubit `q` of a code with parity check tableau `H = parity_checks(code)`
and faults matrix `O = faults_matrix(code)`:

- the flipped detectors are the rows of `H` that anticommute with `P`
  (detector `Dᵢ` ↔ 1-based row `i+1` of `H`, following QuantumClifford's stabilizer row order);
- the flipped logical observables are given by `O * stab_to_gf2(P)` (observable `Lⱼ` ↔
  1-based row `j+1` of `O`, i.e. the first `k` observables are logical-X and the next `k`
  are logical-Z observables). Internally this is read off the columns of `O`: column `q`
  is the image of an `X_q` fault, column `n+q` the image of a `Z_q` fault, and their `xor`
  the image of a `Y_q` fault.

The construction works for any (CSS or non-CSS) stabilizer code for which `parity_checks`
and `faults_matrix` are available, since both the detector and the observable maps are
defined purely through (anti)commutation. The output is deterministic: mechanisms are
emitted qubit by qubit, in `X`, `Y`, `Z` order within a qubit, with sorted target lists.
Faults that flip no detector and no observable (single-qubit Pauli operators that are
elements of the stabilizer group, possible only in rather degenerate codes) carry no
decoding information and are omitted. Faults that flip no detector but do flip an
observable (undetectable logical errors, e.g. any `Z` fault on [`Bitflip3`](@ref)) are kept.

```jldoctest
julia> using QuantumClifford.ECC

julia> dem = detector_error_model(Bitflip3(); px=0.01, pz=0.05);

julia> dem.num_detectors, dem.num_observables, length(dem.errors)
(2, 2, 6)

julia> write_detector_error_model(stdout, dem)
detector D0
detector D1
logical_observable L0
logical_observable L1
error(0.01) D0
error(0.05) L0
error(0.01) D0 D1
error(0.05) L0
error(0.01) D1 L1
error(0.05) L0
```

See the ["ECC example with Pauli Frames"](@ref noisycircuits_pf_ecc_example) page of the documentation
for a complete example, including parsing the exported file with Stim.

See also: [`write_detector_error_model`](@ref), [`faults_matrix`](@ref), [`parity_checks`](@ref)
"""
function detector_error_model(H::Stabilizer; px::Real=0.0, py::Real=0.0, pz::Real=0.0)
    for (name, p) in ((:px, px), (:py, py), (:pz, pz))
        0 <= p <= 1 || throw(ArgumentError("the single-qubit fault probability `$name` must be in [0, 1], got $p"))
    end
    s, n = size(H)
    O = faults_matrix(H)
    num_observables = size(O, 1)
    errors = @NamedTuple{probability::Float64, detectors::Vector{Int}, observables::Vector{Int}}[]
    for q in 1:n
        for (p, has_x, has_z) in ((px, true, false), (py, true, true), (pz, false, true))
            iszero(p) && continue
            fault = has_x ? (has_z ? single_y(n, q) : single_x(n, q)) : single_z(n, q)
            detectors = Int[i-1 for i in 1:s if comm(fault, H, i) == 0x1]
            # equivalent to `O * stab_to_gf2(fault) .% 2` -- see the docstring
            observables = Int[j-1 for j in 1:num_observables if (has_x && O[j, q]) ⊻ (has_z && O[j, n+q])]
            isempty(detectors) && isempty(observables) && continue
            push!(errors, (probability=Float64(p), detectors=detectors, observables=observables))
        end
    end
    return DetectorErrorModel(s, num_observables, errors)
end

detector_error_model(code::AbstractECC; kwargs...) = detector_error_model(parity_checks(code); kwargs...)

"""
$TYPEDSIGNATURES

Write a [`DetectorErrorModel`](@ref) to `io` in Stim's detector error model (`.dem`) text format.

The output starts with explicit `detector D...` and `logical_observable L...` declarations
(so that all detectors and observables are present in the file even if no error mechanism
flips them), followed by one `error(p) D... L...` instruction per error mechanism. The
output is deterministic -- byte-for-byte stable across runs for the same model -- which
makes the exported files easy to test, diff, and cache.

```julia
julia> using QuantumClifford.ECC

julia> dem = detector_error_model(Steane7(); px=1e-3, py=0.0, pz=1e-3);

julia> write_detector_error_model("steane7.dem", dem)
```

The resulting file can be parsed by Stim:

```python
>>> import stim
>>> dem = stim.DetectorErrorModel.from_file("steane7.dem")
>>> dem.num_detectors, dem.num_observables
(6, 2)
```

See also: [`detector_error_model`](@ref)
"""
function write_detector_error_model(io::IO, dem::DetectorErrorModel)
    for d in 0:dem.num_detectors-1
        println(io, "detector D", d)
    end
    for l in 0:dem.num_observables-1
        println(io, "logical_observable L", l)
    end
    for e in dem.errors
        print(io, "error(", e.probability, ')')
        for d in e.detectors
            print(io, " D", d)
        end
        for l in e.observables
            print(io, " L", l)
        end
        println(io)
    end
    return nothing
end

"""
$TYPEDSIGNATURES

Write a [`DetectorErrorModel`](@ref) to the file at `path` in Stim's `.dem` text format.
"""
function write_detector_error_model(path::AbstractString, dem::DetectorErrorModel)
    open(io -> write_detector_error_model(io, dem), path, "w")
    return nothing
end
