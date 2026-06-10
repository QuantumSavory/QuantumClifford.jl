# Importing Stim detector error models (`.dem` files) as Pauli frame circuits.
# Every `error` line becomes an independent `DetectorError` operation and the
# sampling is left to `pftrajectories`. The format is documented at
# https://github.com/quantumlib/Stim/blob/main/doc/file_format_dem_detector_error_model.md

"""
$(TYPEDEF)

A single error mechanism of a detector error model. For each Pauli frame it fires
with probability `p`, and when it does it flips the measurement bits in `bits`.

The entries of `bits` index columns of the [`PauliFrame`](@ref) measurement
buffer. [`read_detector_error_model`](@ref) places the detectors in the first
columns and the logical observables after them, so one `DetectorError` can flip
both kinds of bit.

See also: [`read_detector_error_model`](@ref), [`DeclareMeasurementBits`](@ref).
"""
struct DetectorError{N} <: AbstractOperation
    "The probability that the mechanism fires."
    p::Float64
    "The measurement bit columns flipped when it fires."
    bits::NTuple{N,Int}
end

DetectorError(p, bits::Tuple) = DetectorError{length(bits)}(float(p), bits)
DetectorError(p, detector_bits, logical_bits) = DetectorError(p, (detector_bits..., logical_bits...))

function apply!(frame::PauliFrame, op::DetectorError)
    p = op.p
    @inbounds for f in eachindex(frame)
        if rand() < p
            for b in op.bits
                frame.measurements[f, b] ⊻= true
            end
        end
    end
    return frame
end

# The reference trajectory of a Register is noiseless, so the mechanism never fires there.
apply!(r::Register, ::DetectorError) = r

affectedqubits(::DetectorError) = ()
affectedbits(op::DetectorError) = op.bits

"""
$(TYPEDEF)

A no-op that reserves `nbits` columns in the [`PauliFrame`](@ref) measurement
buffer. [`read_detector_error_model`](@ref) appends one so that detectors or
observables that are declared but never flipped still get an output column.
"""
struct DeclareMeasurementBits <: AbstractOperation
    nbits::Int
end

apply!(frame::PauliFrame, ::DeclareMeasurementBits) = frame
apply!(r::Register, ::DeclareMeasurementBits) = r

affectedqubits(::DeclareMeasurementBits) = ()
affectedbits(op::DeclareMeasurementBits) = (op.nbits,)

"""
$(TYPEDEF)

The circuit obtained from a Stim detector error model. It is a vector of
operations (so it can be passed straight to [`pftrajectories`](@ref)) that also
remembers the number of detectors and observables in the model.

After sampling, `measurements(frame)` is a `trajectories × (num_detectors +
num_observables)` `Bool` matrix. The first `num_detectors` columns are the
detector outcomes `D0, D1, …` and the rest are the observable outcomes
`L0, L1, …`. As always for Pauli frames, these are relative to the noiseless
reference run.

See also: [`read_detector_error_model`](@ref).
"""
struct DetectorErrorModelCircuit <: AbstractVector{AbstractOperation}
    circuit::Vector{AbstractOperation}
    num_detectors::Int
    num_observables::Int
end

Base.size(c::DetectorErrorModelCircuit) = size(c.circuit)
Base.getindex(c::DetectorErrorModelCircuit, i::Int) = c.circuit[i]
Base.IndexStyle(::Type{DetectorErrorModelCircuit}) = IndexLinear()

function Base.show(io::IO, ::MIME"text/plain", c::DetectorErrorModelCircuit)
    nerr = count(op -> op isa DetectorError, c.circuit)
    print(io, "DetectorErrorModelCircuit($(c.num_detectors) detectors, $(c.num_observables) observables, $nerr error mechanisms)")
end

"""
$(TYPEDSIGNATURES)

Import a Stim detector error model (`.dem`) file and return a
[`DetectorErrorModelCircuit`](@ref) that can be sampled with
[`pftrajectories`](@ref).

```julia
using QuantumClifford

circuit = read_detector_error_model("surface_code.dem")
frames = pftrajectories(circuit; trajectories=10_000)
measurements(frames) # trajectories × (num_detectors + num_observables) Bool matrix
```

The supported instructions are `error`, `detector`, `logical_observable`,
`shift_detectors` and (possibly nested) `repeat` blocks, together with comments
and blank lines. The `^` decomposition separators in an `error` line are ignored
for sampling (the mechanism flips every listed target) and coordinates on
`detector`/`shift_detectors` lines are ignored too. Unsupported instructions
raise an `ArgumentError`.

See also: [`detector_error_model_circuit`](@ref), [`parse_detector_error_model`](@ref).
"""
read_detector_error_model(path::AbstractString) = parse_detector_error_model(read(path, String))

"""
$(TYPEDSIGNATURES)

Alias for [`read_detector_error_model`](@ref).
"""
detector_error_model_circuit(path::AbstractString) = read_detector_error_model(path)

"""
$(TYPEDSIGNATURES)

Parse the text of a Stim detector error model into a
[`DetectorErrorModelCircuit`](@ref). See [`read_detector_error_model`](@ref) for
the supported syntax.
"""
function parse_detector_error_model(text::AbstractString)
    lines = _dem_logical_lines(text)
    nodes, _ = _dem_parse_nodes(lines, 1, 0)
    state = _DemState()
    _dem_execute!(nodes, state)
    return _dem_build_circuit(state)
end

# Detector targets are relative to a running offset (advanced by `shift_detectors`);
# observable targets are absolute. We only need the maxima to size the output buffer.
mutable struct _DemState
    detector_offset::Int
    mechanisms::Vector{Tuple{Float64,Vector{Int},Vector{Int}}} # probability, detectors, observables
    max_detector::Int
    max_observable::Int
end

_DemState() = _DemState(0, Tuple{Float64,Vector{Int},Vector{Int}}[], -1, -1)

# Strip comments, drop blanks, and split off the `{` `}` block delimiters (which
# only ever appear in `repeat`) so each lands on its own logical line.
function _dem_logical_lines(text)
    out = String[]
    for raw in split(text, '\n')
        c = findfirst(==('#'), raw)
        code = strip(c === nothing ? raw : raw[1:prevind(raw, c)])
        isempty(code) && continue
        for piece in split(replace(code, '{' => "\n{\n", '}' => "\n}\n"), '\n')
            piece = strip(piece)
            isempty(piece) || push!(out, String(piece))
        end
    end
    return out
end

# Group the logical lines into a tree, turning each `repeat K { … }` into a
# `(K, body)` tuple. `depth` lets us flag unbalanced braces.
function _dem_parse_nodes(lines, i, depth)
    nodes = Any[]
    while i <= length(lines)
        line = lines[i]
        if line == "}"
            depth == 0 && throw(ArgumentError("unmatched '}' in detector error model"))
            return nodes, i + 1
        elseif line == "{"
            throw(ArgumentError("unexpected '{' in detector error model"))
        elseif _dem_head(line) == "repeat"
            m = match(r"^repeat\s+(\d+)$"i, line)
            m === nothing && throw(ArgumentError("malformed repeat instruction: \"$line\""))
            (i < length(lines) && lines[i+1] == "{") || throw(ArgumentError("expected '{' after \"$line\""))
            body, i = _dem_parse_nodes(lines, i + 2, depth + 1)
            push!(nodes, (parse(Int, m.captures[1]), body))
        else
            push!(nodes, line)
            i += 1
        end
    end
    depth == 0 || throw(ArgumentError("unmatched '{' in detector error model"))
    return nodes, i
end

function _dem_head(line)
    m = match(r"^[A-Za-z_]+", line)
    m === nothing && throw(ArgumentError("cannot parse detector error model line: \"$line\""))
    return lowercase(m.match)
end

function _dem_execute!(nodes, state)
    for node in nodes
        if node isa Tuple
            reps, body = node
            for _ in 1:reps
                _dem_execute!(body, state)
            end
        else
            _dem_exec_instruction!(node, state)
        end
    end
    return state
end

function _dem_exec_instruction!(line, state)
    head = _dem_head(line)
    args, targets = _dem_args_targets(line)
    if head == "error"
        length(args) == 1 || throw(ArgumentError("`error` needs exactly one probability argument in line: \"$line\""))
        p = _dem_parse_float(args[1], line)
        0 <= p <= 1 || throw(ArgumentError("error probability must be in [0, 1] but is $p in line: \"$line\""))
        dets, obs = Int[], Int[]
        for t in targets
            _dem_collect_target!(t, state, dets, obs, line)
        end
        isempty(dets) || (state.max_detector = max(state.max_detector, maximum(dets)))
        isempty(obs) || (state.max_observable = max(state.max_observable, maximum(obs)))
        (isempty(dets) && isempty(obs)) || push!(state.mechanisms, (p, dets, obs))
    elseif head == "detector"
        for t in targets
            state.max_detector = max(state.max_detector, _dem_detector_index(t, state, line))
        end
    elseif head == "logical_observable"
        for t in targets
            state.max_observable = max(state.max_observable, _dem_observable_index(t, line))
        end
    elseif head == "shift_detectors"
        length(targets) == 1 || throw(ArgumentError("`shift_detectors` needs exactly one integer target in line: \"$line\""))
        state.detector_offset += _dem_parse_int(targets[1], line)
    else
        throw(ArgumentError("unsupported detector error model instruction \"$head\" in line: \"$line\""))
    end
    return state
end

# Separate the parenthesized arguments from the space-separated targets, skipping
# the optional `[tag]`. Coordinate arguments are returned unparsed and unused.
function _dem_args_targets(line)
    name = match(r"^[A-Za-z_]+", line)
    rest = strip(line[nextind(line, name.offset + length(name.match) - 1):end])
    if startswith(rest, '[')
        b = findfirst(==(']'), rest)
        b === nothing && throw(ArgumentError("unterminated tag in line: \"$line\""))
        rest = strip(rest[nextind(rest, b):end])
    end
    args = SubString[]
    if startswith(rest, '(')
        b = findfirst(==(')'), rest)
        b === nothing && throw(ArgumentError("unterminated '(' in line: \"$line\""))
        args = [strip(a) for a in split(rest[2:prevind(rest, b)], ',') if !isempty(strip(a))]
        rest = strip(rest[nextind(rest, b):end])
    end
    targets = isempty(rest) ? String[] : split(rest)
    return args, targets
end

function _dem_collect_target!(t, state, dets, obs, line)
    t == "^" && return # decomposition hint, irrelevant for sampling
    if first(t) in ('D', 'd')
        push!(dets, state.detector_offset + _dem_parse_int(t[2:end], line))
    elseif first(t) in ('L', 'l')
        push!(obs, _dem_parse_int(t[2:end], line))
    else
        throw(ArgumentError("unsupported target \"$t\" in detector error model line: \"$line\""))
    end
    return nothing
end

function _dem_detector_index(t, state, line)
    first(t) in ('D', 'd') || throw(ArgumentError("expected a detector target (D#) but found \"$t\" in line: \"$line\""))
    return state.detector_offset + _dem_parse_int(t[2:end], line)
end

function _dem_observable_index(t, line)
    first(t) in ('L', 'l') || throw(ArgumentError("expected an observable target (L#) but found \"$t\" in line: \"$line\""))
    return _dem_parse_int(t[2:end], line)
end

function _dem_parse_int(s, line)
    v = tryparse(Int, strip(String(s)))
    v === nothing && throw(ArgumentError("expected an integer but found \"$s\" in detector error model line: \"$line\""))
    return v
end

function _dem_parse_float(s, line)
    v = tryparse(Float64, strip(String(s)))
    v === nothing && throw(ArgumentError("expected a number but found \"$s\" in detector error model line: \"$line\""))
    return v
end

# Detectors occupy columns 1:D, observables D+1:D+L.
function _dem_build_circuit(state)
    D = state.max_detector + 1
    L = state.max_observable + 1
    D + L == 0 && throw(ArgumentError("the detector error model declares no detectors or logical observables"))
    circuit = AbstractOperation[]
    for (p, dets, obs) in state.mechanisms
        push!(circuit, DetectorError(p, [d + 1 for d in dets], [D + o + 1 for o in obs]))
    end
    push!(circuit, DeclareMeasurementBits(D + L))
    return DetectorErrorModelCircuit(circuit, D, L)
end
