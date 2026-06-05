##############################
# Stim detector error model (.dem) import
#
# This file implements importing of Stim "detector error model" (`.dem`) files
# as QuantumClifford circuits that can be sampled with `pftrajectories`.
#
# A detector error model is a list of independent error mechanisms. Each
# mechanism has a probability `p` and a set of *symptoms* (detectors) and
# *frame changes* (logical observables) it flips when it fires. Sampling such a
# model is exactly a Pauli-frame style simulation: for every trajectory and
# every mechanism, flip a biased coin and, if it comes up heads, XOR the listed
# detector and logical bits into the trajectory's measurement record.
#
# The detector and logical outcomes are stored in the `measurements` matrix of
# a `PauliFrame`. The column layout is:
#   * columns `1:num_detectors`                     -> detectors `D0, D1, ...`
#   * columns `num_detectors .+ (1:num_logicals)`   -> observables `L0, L1, ...`
#
# Reference for the file format:
# https://github.com/quantumlib/Stim/blob/main/doc/file_format_dem_detector_error_model.md
##############################

"""
$(TYPEDEF)

A single independent error mechanism of a Stim detector error model, lowered to
a [`PauliFrame`](@ref) operation.

When applied to a [`PauliFrame`](@ref), an independent Bernoulli event with
probability `p` is sampled for every trajectory. Whenever the event occurs, all
the measurement-storage columns listed in `bits` are flipped (XOR-ed with `1`)
simultaneously. The `bits` are absolute column indices into the measurement
storage of a [`PauliFrame`](@ref); the mapping from detectors/logical
observables to columns is established by [`read_detector_error_model`](@ref).

This is the operation that the Stim `error(p) D... L...` instruction is lowered
to. See also: [`read_detector_error_model`](@ref).
"""
struct DetectorError <: AbstractOperation
    "the probability that this error mechanism fires in a given trajectory"
    p::Float64
    "the measurement-storage columns (detectors and observables) flipped together when the mechanism fires"
    bits::Vector{Int}
end

"""
$(TYPEDEF)

A no-op [`PauliFrame`](@ref) operation that only serves to reserve
measurement-storage columns.

Detectors and logical observables can be declared in a detector error model
without ever being touched by an error mechanism. Such columns still need to be
allocated in the output of [`pftrajectories`](@ref). This operation carries the
full list of columns of an imported model so that the measurement matrix is
always sized to `num_detectors + num_logicals`, regardless of which columns are
actually flipped by some [`DetectorError`](@ref).
"""
struct DetectorErrorModelColumns <: AbstractOperation
    "all measurement-storage columns reserved by the imported detector error model"
    bits::Vector{Int}
end

affectedqubits(::DetectorError) = ()
affectedbits(op::DetectorError) = op.bits

affectedqubits(::DetectorErrorModelColumns) = ()
affectedbits(op::DetectorErrorModelColumns) = op.bits

function apply!(frame::PauliFrame, op::DetectorError)
    p = op.p
    bits = op.bits
    isempty(bits) && return frame
    meas = frame.measurements
    @inbounds for f in eachindex(frame)
        if rand() < p
            for b in bits
                meas[f, b] ⊻= true
            end
        end
    end
    return frame
end

apply!(frame::PauliFrame, ::DetectorErrorModelColumns) = frame

"""
$(TYPEDEF)

A circuit-like container holding the operations of an imported Stim detector
error model together with the number of detectors and logical observables it
declares.

It behaves like a vector of [`AbstractOperation`](@ref) (so it can be passed
directly to [`pftrajectories`](@ref)) while also recording `num_detectors` and
`num_logicals`, which describe the column layout of the resulting measurement
matrix:

  * columns `1:num_detectors` hold the detector outcomes `D0, D1, ...`
  * columns `num_detectors .+ (1:num_logicals)` hold the observable outcomes `L0, L1, ...`

See also: [`read_detector_error_model`](@ref).
"""
struct DetectorErrorModelCircuit <: AbstractVector{AbstractOperation}
    ops::Vector{AbstractOperation}
    num_detectors::Int
    num_logicals::Int
end

Base.size(c::DetectorErrorModelCircuit) = size(c.ops)
Base.getindex(c::DetectorErrorModelCircuit, i::Int) = c.ops[i]
Base.IndexStyle(::Type{DetectorErrorModelCircuit}) = IndexLinear()

function Base.show(io::IO, ::MIME"text/plain", c::DetectorErrorModelCircuit)
    nmech = count(o->isa(o, DetectorError), c.ops)
    print(io, "DetectorErrorModelCircuit with $(c.num_detectors) detector(s), ",
              "$(c.num_logicals) logical observable(s), and $(nmech) error mechanism(s)")
end

##############################
# Parsing
##############################

# Mutable state tracked while interpreting a `.dem` file.
mutable struct _DEMState
    detector_offset::Int                                 # current relative detector index offset
    mechanisms::Vector{Tuple{Float64,Vector{Int},Vector{Int}}} # (p, abs detector indices, logical indices), 0-based
    max_detector::Int                                    # largest absolute detector index seen (-1 if none)
    max_logical::Int                                     # largest logical index seen (-1 if none)
end
_DEMState() = _DEMState(0, Tuple{Float64,Vector{Int},Vector{Int}}[], -1, -1)

# Remove a trailing `# ...` comment, while not treating a `#` inside a `[tag]` as a comment.
function _strip_comment(line::AbstractString)
    inbracket = false
    for (i, c) in pairs(line)
        if c == '['
            inbracket = true
        elseif c == ']'
            inbracket = false
        elseif c == '#' && !inbracket
            return SubString(line, 1, prevind(line, i))
        end
    end
    return SubString(line, 1)
end

# Parse a single instruction line (already stripped of comments and surrounding whitespace)
# into `(name, args, targets, isblock)`.
function _parse_instruction(line::AbstractString)
    isblock = false
    if endswith(line, "{")
        isblock = true
        line = strip(SubString(line, 1, prevind(line, lastindex(line))))
    end
    m = match(r"^([A-Za-z][A-Za-z0-9_]*)", line)
    m === nothing && throw(ArgumentError(lazy"Could not parse a detector-error-model instruction name in line: `$line`"))
    name = lowercase(m.match)
    rest = strip(SubString(line, ncodeunits(m.match) + 1))
    # optional tag, e.g. `error[mytag](...)` -- parsed and discarded
    if startswith(rest, "[")
        idx = findfirst(']', rest)
        idx === nothing && throw(ArgumentError(lazy"Unterminated instruction tag `[` in line: `$line`"))
        rest = strip(SubString(rest, idx + 1))
    end
    args = Float64[]
    if startswith(rest, "(")
        idx = findfirst(')', rest)
        idx === nothing && throw(ArgumentError(lazy"Unterminated argument list `(` in line: `$line`"))
        argstr = SubString(rest, 2, prevind(rest, idx))
        for a in split(argstr, ',')
            a = strip(a)
            isempty(a) && continue
            v = tryparse(Float64, a)
            v === nothing && throw(ArgumentError(lazy"Could not parse numeric argument `$a` in line: `$line`"))
            push!(args, v)
        end
        rest = strip(SubString(rest, idx + 1))
    end
    targets = isempty(rest) ? SubString{String}[] : split(rest)
    return name, args, targets, isblock
end

# Parse a single target token into `(kind, index)` where kind is one of
# `:D` (detector), `:L` (logical), `:sep` (decomposition separator), `:num` (numeric).
function _parse_target(t::AbstractString, line)
    if (startswith(t, "D") || startswith(t, "d")) && length(t) > 1
        idx = tryparse(Int, SubString(t, 2))
        (idx === nothing || idx < 0) && throw(ArgumentError(lazy"Invalid detector target `$t` in line: `$line`"))
        return :D, idx
    elseif (startswith(t, "L") || startswith(t, "l")) && length(t) > 1
        idx = tryparse(Int, SubString(t, 2))
        (idx === nothing || idx < 0) && throw(ArgumentError(lazy"Invalid logical-observable target `$t` in line: `$line`"))
        return :L, idx
    elseif t == "^"
        return :sep, 0
    else
        idx = tryparse(Int, t)
        idx === nothing && throw(ArgumentError(lazy"Invalid target `$t` in line: `$line`"))
        return :num, idx
    end
end

# Toggle membership of `x` in `v` (parity tracking): if present, remove it; otherwise add it.
function _toggle!(v::Vector{Int}, x::Int)
    i = findfirst(==(x), v)
    if i === nothing
        push!(v, x)
    else
        deleteat!(v, i)
    end
    return v
end

function _apply_instruction!(state::_DEMState, name, args, targets, line)
    if name == "error"
        length(args) == 1 || throw(ArgumentError(lazy"`error` requires exactly one probability argument, got $(length(args)), in line: `$line`"))
        p = args[1]
        (0.0 <= p <= 1.0) || throw(ArgumentError(lazy"`error` probability must lie in [0,1], got $p, in line: `$line`"))
        dets = Int[]
        logs = Int[]
        for t in targets
            kind, idx = _parse_target(t, line)
            if kind === :sep # decomposition hint -- only the overall (parity) symptom matters for sampling
                continue
            elseif kind === :D
                abs = idx + state.detector_offset
                _toggle!(dets, abs)
                state.max_detector = max(state.max_detector, abs)
            elseif kind === :L
                _toggle!(logs, idx)
                state.max_logical = max(state.max_logical, idx)
            else
                throw(ArgumentError(lazy"`error` targets must be detectors (D#), logical observables (L#), or separators (^); got `$t` in line: `$line`"))
            end
        end
        push!(state.mechanisms, (p, dets, logs))
    elseif name == "detector"
        isempty(targets) && throw(ArgumentError(lazy"`detector` instruction requires at least one detector target in line: `$line`"))
        for t in targets
            kind, idx = _parse_target(t, line)
            kind === :D || throw(ArgumentError(lazy"`detector` targets must be detectors (D#); got `$t` in line: `$line`"))
            abs = idx + state.detector_offset
            state.max_detector = max(state.max_detector, abs)
        end
    elseif name == "logical_observable"
        isempty(targets) && throw(ArgumentError(lazy"`logical_observable` instruction requires at least one logical-observable target in line: `$line`"))
        for t in targets
            kind, idx = _parse_target(t, line)
            kind === :L || throw(ArgumentError(lazy"`logical_observable` targets must be logical observables (L#); got `$t` in line: `$line`"))
            state.max_logical = max(state.max_logical, idx)
        end
    elseif name == "shift_detectors"
        length(targets) <= 1 || throw(ArgumentError(lazy"`shift_detectors` takes at most one numeric target, got $(length(targets)), in line: `$line`"))
        if length(targets) == 1
            kind, inc = _parse_target(targets[1], line)
            (kind === :num && inc >= 0) || throw(ArgumentError(lazy"`shift_detectors` numeric target must be a non-negative integer; got `$(targets[1])` in line: `$line`"))
            state.detector_offset += inc
        end
    else
        throw(ArgumentError(lazy"Unsupported detector-error-model instruction `$name` in line: `$line`"))
    end
    return state
end

# Find the index just past the `}` matching a block that starts at `start`.
function _find_block_end(lines, start::Int)
    depth = 1
    i = start
    n = length(lines)
    while i <= n
        line = strip(_strip_comment(lines[i]))
        if line == "}"
            depth -= 1
            depth == 0 && return i + 1
        elseif endswith(line, "{")
            depth += 1
        end
        i += 1
    end
    throw(ArgumentError("Unterminated `repeat` block: missing `}`."))
end

# Interpret `lines` starting at `start`. Returns the index just past where this
# block stops: at end of input for the top level, or just past the matching `}`.
function _run_block!(state::_DEMState, lines, start::Int, toplevel::Bool)
    i = start
    n = length(lines)
    while i <= n
        line = strip(_strip_comment(lines[i]))
        if isempty(line)
            i += 1
            continue
        end
        if line == "}"
            toplevel && throw(ArgumentError("Unexpected `}` with no matching block opener."))
            return i + 1
        end
        name, args, targets, isblock = _parse_instruction(line)
        if isblock
            name == "repeat" || throw(ArgumentError(lazy"Unsupported block instruction `$name` in line: `$line`"))
            length(targets) == 1 || throw(ArgumentError(lazy"`repeat` requires exactly one numeric repetition count in line: `$line`"))
            kind, k = _parse_target(targets[1], line)
            (kind === :num && k >= 0) || throw(ArgumentError(lazy"Invalid `repeat` count `$(targets[1])` in line: `$line`"))
            body_start = i + 1
            block_end = body_start
            for _ in 1:k
                block_end = _run_block!(state, lines, body_start, false)
            end
            k == 0 && (block_end = _find_block_end(lines, body_start))
            i = block_end
        else
            _apply_instruction!(state, name, args, targets, line)
            i += 1
        end
    end
    toplevel || throw(ArgumentError("Unterminated `repeat` block: missing `}`."))
    return i
end

function _build_circuit(state::_DEMState)
    D = state.max_detector + 1 # number of detector columns
    L = state.max_logical + 1  # number of logical-observable columns
    ncols = D + L
    ops = AbstractOperation[]
    if ncols > 0
        push!(ops, DetectorErrorModelColumns(collect(1:ncols)))
    end
    for (p, dets, logs) in state.mechanisms
        bits = Int[]
        for d in dets
            push!(bits, d + 1)
        end
        for l in logs
            push!(bits, D + l + 1)
        end
        push!(ops, DetectorError(p, bits))
    end
    return DetectorErrorModelCircuit(ops, D, L)
end

_parse_detector_error_model(lines::AbstractVector) = begin
    state = _DEMState()
    _run_block!(state, lines, 1, true)
    _build_circuit(state)
end

"""
$(TYPEDSIGNATURES)

Import a Stim detector error model (`.dem`) file and return a circuit-like
object (a [`DetectorErrorModelCircuit`](@ref)) that can be sampled directly with
[`pftrajectories`](@ref).

A detector error model is a list of independent error mechanisms. Each
`error(p) D... L...` instruction becomes one [`DetectorError`](@ref) operation
which, with probability `p`, flips the listed detector and logical-observable
bits in a trajectory's measurement record. Sampling is performed entirely by the
existing [`pftrajectories`](@ref) backend.

The supported subset of the [`.dem` format](https://github.com/quantumlib/Stim/blob/main/doc/file_format_dem_detector_error_model.md)
includes the `error`, `detector`, `logical_observable`, and `shift_detectors`
instructions, `repeat K { ... }` blocks (including nesting), comments, and blank
lines. The `^` decomposition separators inside `error` instructions are
accepted; since they only suggest a decomposition, the importer keeps the
overall (parity) set of symptoms and frame changes. Detector and observable
coordinate annotations are parsed but not stored.

The returned object indexes the measurement matrix of the resulting
[`PauliFrame`](@ref) as follows:

  * columns `1:num_detectors` hold detector outcomes `D0, D1, ...`
  * columns `num_detectors .+ (1:num_logicals)` hold observable outcomes `L0, L1, ...`

where `num_detectors` and `num_logicals` are fields of the returned
[`DetectorErrorModelCircuit`](@ref).

# Example

```julia
circuit = read_detector_error_model("surface_code.dem")
frames  = pftrajectories(circuit; trajectories=10_000)
m = measurements(frames)                                  # trajectories × (num_detectors+num_logicals)
detector_samples = m[:, 1:circuit.num_detectors]
logical_samples  = m[:, circuit.num_detectors .+ (1:circuit.num_logicals)]
```

See also: [`pftrajectories`](@ref), [`DetectorError`](@ref).
"""
function read_detector_error_model(path::AbstractString)
    lines = readlines(path)
    return _parse_detector_error_model(lines)
end

"""
$(TYPEDSIGNATURES)

Import a Stim detector error model directly from an `IO` stream. See
[`read_detector_error_model(::AbstractString)`](@ref) for details.
"""
function read_detector_error_model(io::IO)
    lines = readlines(io)
    return _parse_detector_error_model(lines)
end
