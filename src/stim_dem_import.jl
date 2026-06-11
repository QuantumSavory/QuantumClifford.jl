# Importing Stim detector error model (.dem) files as circuits of operations
# that can be sampled with `pftrajectories`.
#
# See `read_detector_error_model`, `parse_detector_error_model`,
# `DetectorError`, `DemDeclaration`, and `DetectorErrorModelCircuit`.
#
# The format is specified in the Stim repository:
# https://github.com/quantumlib/Stim/blob/main/doc/file_format_dem_detector_error_model.md

##############################
# Operations
##############################

"""
$(TYPEDEF)

An independent classical error mechanism from a detector error model.

With probability `p` (independently for every frame/trajectory), the error occurs
and flips (XORs) the measurement bits listed in `detector_bits` and `logical_bits`.
The two lists are kept separate only for bookkeeping purposes (they are treated
identically during sampling): when a detector error model is imported with
[`read_detector_error_model`](@ref), detectors are assigned the first block of
measurement-bit columns and logical observables the following block.

This operation acts only on the classical measurement record -- it affects no qubits.
It corresponds to the `error(p) D... L...` instruction of Stim's `.dem` format.

```jldoctest
julia> circuit = [DemDeclaration(2, 1), DetectorError(1.0, [1, 2], [3])];

julia> frames = pftrajectories(circuit; trajectories=4, threads=false);

julia> measurements(frames)
4×3 Matrix{Bool}:
 1  1  1
 1  1  1
 1  1  1
 1  1  1
```

See also: [`read_detector_error_model`](@ref), [`DemDeclaration`](@ref)
"""
struct DetectorError <: AbstractOperation
    "The probability with which this error mechanism occurs (independently per trajectory)."
    p::Float64
    "Measurement-bit (column) indices, 1-based, that store detector outcomes flipped by this error."
    detector_bits::Vector{Int}
    "Measurement-bit (column) indices, 1-based, that store logical-observable outcomes flipped by this error."
    logical_bits::Vector{Int}
    function DetectorError(p, detector_bits, logical_bits)
        0 <= p <= 1 || throw(ArgumentError(lazy"`DetectorError` requires a probability (0 to 1) but got $(p)."))
        d = collect(Int, detector_bits)
        l = collect(Int, logical_bits)
        (all(>(0), d) && all(>(0), l)) || throw(ArgumentError("`DetectorError` measurement-bit indices have to be larger than zero. Ensure indexing always starts from 1."))
        new(Float64(p), d, l)
    end
end

"""
$(TYPEDEF)

A no-op operation declaring how many detectors and logical observables an imported
detector error model contains.

Its only effect is through its `affectedbits` method, which reports the full range
`1:(n_detectors+n_logicals)` of measurement-bit columns, so that
[`pftrajectories`](@ref) allocates measurement storage for *all* declared detectors
and logical observables -- including ones that no [`DetectorError`](@ref) ever flips
(e.g. detectors only declared via a `detector` instruction).
[`read_detector_error_model`](@ref) places one `DemDeclaration` at the start of the
returned circuit.

See also: [`read_detector_error_model`](@ref), [`DetectorError`](@ref)
"""
struct DemDeclaration <: AbstractOperation
    "Number of detectors (occupying measurement-bit columns `1:n_detectors`)."
    n_detectors::Int
    "Number of logical observables (occupying measurement-bit columns `n_detectors+1:n_detectors+n_logicals`)."
    n_logicals::Int
    function DemDeclaration(n_detectors, n_logicals)
        (n_detectors >= 0 && n_logicals >= 0) || throw(ArgumentError("`DemDeclaration` counts must be non-negative."))
        new(n_detectors, n_logicals)
    end
end

Base.:(==)(l::DetectorError, r::DetectorError) = l.p == r.p && l.detector_bits == r.detector_bits && l.logical_bits == r.logical_bits
Base.hash(d::DetectorError, h::UInt) = hash((DetectorError, d.p, d.detector_bits, d.logical_bits), h)

affectedqubits(::DetectorError) = ()
affectedqubits(::DemDeclaration) = ()
affectedbits(d::DetectorError) = vcat(d.detector_bits, d.logical_bits)
affectedbits(d::DemDeclaration) = 1:(d.n_detectors + d.n_logicals)

function apply!(frame::PauliFrame, op::DetectorError)
    p = op.p
    bits = frame.measurements
    for f in eachindex(frame)
        if rand() < p
            for b in op.detector_bits
                bits[f, b] ⊻= true
            end
            for b in op.logical_bits
                bits[f, b] ⊻= true
            end
        end
    end
    return frame
end

apply!(frame::PauliFrame, ::DemDeclaration) = frame

function apply!(r::Register, op::DetectorError)
    if rand() < op.p
        bits = bitview(r)
        for b in op.detector_bits
            bits[b] ⊻= true
        end
        for b in op.logical_bits
            bits[b] ⊻= true
        end
    end
    return r
end

apply!(r::Register, ::DemDeclaration) = r

##############################
# The circuit-like container
##############################

"""
$(TYPEDEF)

A circuit imported from a Stim detector error model (`.dem`) file -- a plain vector
of operations ([`DemDeclaration`](@ref) followed by [`DetectorError`](@ref)s)
together with the number of detectors and logical observables.

It is an `AbstractVector{AbstractOperation}`, so it is accepted directly by the
existing [`pftrajectories`](@ref) and [`mctrajectories`](@ref) methods.

The measurement record sampled from this circuit has one row per trajectory and
`n_detectors + n_logicals` columns: detector `Dk` of the `.dem` file is stored in
column `k+1` and logical observable `Lk` in column `n_detectors + k + 1`.
Use [`detectorview`](@ref) and [`observableview`](@ref) to slice the two blocks apart.

See also: [`read_detector_error_model`](@ref), [`parse_detector_error_model`](@ref)
"""
struct DetectorErrorModelCircuit <: AbstractVector{AbstractOperation}
    "The operations making up the circuit."
    ops::Vector{AbstractOperation}
    "Number of detectors (first `n_detectors` measurement-bit columns)."
    n_detectors::Int
    "Number of logical observables (next `n_logicals` measurement-bit columns)."
    n_logicals::Int
end

Base.size(c::DetectorErrorModelCircuit) = size(c.ops)
Base.getindex(c::DetectorErrorModelCircuit, i::Int) = c.ops[i]
Base.IndexStyle(::Type{DetectorErrorModelCircuit}) = IndexLinear()
Base.:(==)(l::DetectorErrorModelCircuit, r::DetectorErrorModelCircuit) = l.n_detectors == r.n_detectors && l.n_logicals == r.n_logicals && l.ops == r.ops

"""
$(TYPEDSIGNATURES)

A view of the detector columns of the measurement record sampled from an imported
detector error model circuit: a `trajectories × n_detectors` `Bool` array in which
entry `[i, k+1]` says whether detector `Dk` was flipped in trajectory `i`.

```jldoctest
julia> circuit = parse_detector_error_model("error(1) D1 L0");

julia> frames = pftrajectories(circuit; trajectories=2, threads=false);

julia> collect(detectorview(circuit, frames))
2×2 Matrix{Bool}:
 0  1
 0  1
```

See also: [`observableview`](@ref), [`measurements`](@ref)
"""
detectorview(c::DetectorErrorModelCircuit, m::AbstractMatrix) = view(m, :, 1:c.n_detectors)
detectorview(c::DetectorErrorModelCircuit, f::PauliFrame) = detectorview(c, measurements(f))

"""
$(TYPEDSIGNATURES)

A view of the logical-observable columns of the measurement record sampled from an
imported detector error model circuit: a `trajectories × n_logicals` `Bool` array
in which entry `[i, k+1]` says whether logical observable `Lk` was flipped in
trajectory `i`.

```jldoctest
julia> circuit = parse_detector_error_model("error(1) D1 L0");

julia> frames = pftrajectories(circuit; trajectories=2, threads=false);

julia> collect(observableview(circuit, frames))
2×1 Matrix{Bool}:
 1
 1
```

See also: [`detectorview`](@ref), [`measurements`](@ref)
"""
observableview(c::DetectorErrorModelCircuit, m::AbstractMatrix) = view(m, :, c.n_detectors+1:c.n_detectors+c.n_logicals)
observableview(c::DetectorErrorModelCircuit, f::PauliFrame) = observableview(c, measurements(f))

##############################
# Parsing
##############################

# Internal parsed-instruction representation (a small tree; `repeat` blocks are
# kept as nodes and only unrolled during interpretation).
abstract type AbstractDemInstruction end

struct DemErrorInstruction <: AbstractDemInstruction
    p::Float64
    detectors::Vector{Int} # relative detector ids, in file order, possibly repeating
    logicals::Vector{Int}  # logical observable ids, in file order, possibly repeating
end
struct DemDetectorInstruction <: AbstractDemInstruction
    ids::Vector{Int}       # relative detector ids
end
struct DemLogicalInstruction <: AbstractDemInstruction
    ids::Vector{Int}       # logical observable ids
end
struct DemShiftInstruction <: AbstractDemInstruction
    shift::Int
end
struct DemRepeatInstruction <: AbstractDemInstruction
    count::Int
    body::Vector{AbstractDemInstruction}
end

_dem_err(lineno, msg) = throw(ArgumentError("Detector error model parsing failed on line $(lineno): $(msg)"))

# Remove a `# comment` from a line, being careful that `#` may legally appear
# inside an instruction tag (`error[my#tag](0.1) D0`); `]` cannot appear inside
# a tag (it is escaped as `\\C`), so scanning for the first `]` is correct.
function _dem_strip_comment(line::AbstractString)
    intag = false
    for j in eachindex(line)
        c = line[j]
        if intag
            c == ']' && (intag = false)
        elseif c == '['
            intag = true
        elseif c == '#'
            return line[1:prevind(line, j)]
        end
    end
    return line
end

function _dem_parse_args(lineno, s::AbstractString) # the inside of `(...)`
    pieces = split(s, ',')
    args = Vector{Float64}(undef, length(pieces))
    for (i, piece) in enumerate(pieces)
        v = tryparse(Float64, strip(piece))
        v === nothing && _dem_err(lineno, "expected a number as an instruction argument but got '$(strip(piece))'.")
        args[i] = v
    end
    return args
end

function _dem_parse_uint(lineno, s::AbstractString, what)
    v = isempty(s) ? nothing : tryparse(Int, s)
    (v === nothing || !all(isdigit, s)) && _dem_err(lineno, "expected a non-negative integer as $(what) but got '$(s)'.")
    return v
end

# Parses instructions from `lines` starting at index `i` until end-of-file
# (`inblock=false`) or until a matching `}` (`inblock=true`).
# Returns `(instructions, index_of_next_unread_line)`.
function _dem_parse_block(lines, i, inblock)
    instructions = AbstractDemInstruction[]
    while i <= length(lines)
        lineno = i
        line = strip(_dem_strip_comment(lines[i]))
        i += 1
        isempty(line) && continue
        if line == "}"
            inblock && return instructions, i
            _dem_err(lineno, "uninitiated block: got a '}' without a '{'.")
        end
        # block start? (the `{` is always on the same line as its instruction)
        isblockstart = endswith(line, '{')
        if isblockstart
            line = strip(line[1:prevind(line, lastindex(line))])
        end
        # instruction name (case insensitive)
        m = match(r"^[a-zA-Z][a-zA-Z0-9_]*", line)
        m === nothing && _dem_err(lineno, "expected an instruction name but got '$(line)'.")
        name = lowercase(m.match)
        rest = SubString(line, m.offset + ncodeunits(m.match)) # the name is ASCII, so byte offsets are safe
        # optional [tag] -- parsed and discarded
        if startswith(rest, '[')
            closing = findfirst(']', rest)
            closing === nothing && _dem_err(lineno, "unterminated instruction tag (missing ']').")
            rest = rest[nextind(rest, closing):end]
        end
        # optional (arguments)
        args = Float64[]
        hasargs = false
        if startswith(rest, '(')
            closing = findfirst(')', rest)
            closing === nothing && _dem_err(lineno, "unterminated argument list (missing ')').")
            hasargs = true
            args = _dem_parse_args(lineno, rest[nextind(rest, 1):prevind(rest, closing)])
            rest = rest[nextind(rest, closing):end]
        end
        targets = split(rest)
        if name == "repeat"
            isblockstart || _dem_err(lineno, "missing '{' at start of repeat block.")
            hasargs && _dem_err(lineno, "'repeat' does not take parenthesized arguments.")
            length(targets) == 1 || _dem_err(lineno, "'repeat' takes exactly 1 numeric target (the repetition count) but got $(length(targets)) targets.")
            count = _dem_parse_uint(lineno, targets[1], "the repetition count")
            body, i = _dem_parse_block(lines, i, true)
            push!(instructions, DemRepeatInstruction(count, body))
            continue
        end
        isblockstart && _dem_err(lineno, "unexpected '{' (only 'repeat' instructions start blocks).")
        if name == "error"
            length(args) == 1 || _dem_err(lineno, "'error' takes exactly 1 argument (a probability) but got $(length(args)) arguments.")
            p = args[1]
            0 <= p <= 1 || _dem_err(lineno, "'error' argument must be a probability (0 to 1) but got $(p).")
            detectors = Int[]
            logicals = Int[]
            lastsep = true # also catches a separator as the first target
            for t in targets
                if t == "^"
                    lastsep && _dem_err(lineno, "'error' separators (^) must separate non-empty groups of targets.")
                    lastsep = true
                    continue
                end
                lastsep = false
                if startswith(t, 'D') || startswith(t, 'd')
                    push!(detectors, _dem_parse_uint(lineno, t[2:end], "a detector id"))
                elseif startswith(t, 'L') || startswith(t, 'l')
                    push!(logicals, _dem_parse_uint(lineno, t[2:end], "a logical observable id"))
                else
                    _dem_err(lineno, "unrecognized 'error' target '$(t)' (expected D#, L#, or ^).")
                end
            end
            lastsep && !isempty(targets) && _dem_err(lineno, "'error' separators (^) must separate non-empty groups of targets.")
            push!(instructions, DemErrorInstruction(p, detectors, logicals))
        elseif name == "detector"
            # arguments are coordinates -- parsed and discarded
            isempty(targets) && _dem_err(lineno, "'detector' takes at least 1 target (D#).")
            ids = Int[]
            for t in targets
                (startswith(t, 'D') || startswith(t, 'd')) || _dem_err(lineno, "'detector' takes relative detector targets (D#) but got '$(t)'.")
                push!(ids, _dem_parse_uint(lineno, t[2:end], "a detector id"))
            end
            push!(instructions, DemDetectorInstruction(ids))
        elseif name == "logical_observable"
            hasargs && _dem_err(lineno, "'logical_observable' takes 0 arguments but got $(length(args)) arguments.")
            isempty(targets) && _dem_err(lineno, "'logical_observable' takes at least 1 target (L#).")
            ids = Int[]
            for t in targets
                (startswith(t, 'L') || startswith(t, 'l')) || _dem_err(lineno, "'logical_observable' takes logical observable targets (L#) but got '$(t)'.")
                push!(ids, _dem_parse_uint(lineno, t[2:end], "a logical observable id"))
            end
            push!(instructions, DemLogicalInstruction(ids))
        elseif name == "shift_detectors"
            # arguments are coordinate shifts -- parsed and discarded
            length(targets) == 1 || _dem_err(lineno, "'shift_detectors' takes exactly 1 numeric target but got $(length(targets)) targets.")
            push!(instructions, DemShiftInstruction(_dem_parse_uint(lineno, targets[1], "the detector index shift")))
        else
            _dem_err(lineno, "unrecognized instruction name '$(name)'. Supported instructions: error, detector, logical_observable, shift_detectors, repeat.")
        end
    end
    inblock && throw(ArgumentError("Detector error model parsing failed: unterminated block -- got a '{' without an eventual '}'."))
    return instructions, i
end

##############################
# Interpretation (unrolling)
##############################

mutable struct DemInterpreterState
    offset::Int  # current detector index offset, accumulated from `shift_detectors`
    maxdet::Int  # largest absolute detector index mentioned so far (-1 if none)
    maxobs::Int  # largest logical observable index mentioned so far (-1 if none)
    mechanisms::Vector{Tuple{Float64,Vector{Int},Vector{Int}}} # (p, absolute detector ids, observable ids)
end

# Keep only the ids that appear an odd number of times (Stim semantics: targets
# repeated within one `error` instruction cancel out), sorted.
function _odd_parity(ids)
    counts = Dict{Int,Int}()
    for k in ids
        counts[k] = get(counts, k, 0) + 1
    end
    return sort!([k for (k, c) in counts if isodd(c)])
end

function _dem_interpret!(st::DemInterpreterState, instructions)
    for instruction in instructions
        if instruction isa DemShiftInstruction
            st.offset += instruction.shift
        elseif instruction isa DemDetectorInstruction
            for k in instruction.ids
                st.maxdet = max(st.maxdet, st.offset + k)
            end
        elseif instruction isa DemLogicalInstruction
            for k in instruction.ids
                st.maxobs = max(st.maxobs, k)
            end
        elseif instruction isa DemErrorInstruction
            detectors = [st.offset + k for k in instruction.detectors]
            # All mentioned ids count towards the model size (matching Stim's
            # `count_detectors`/`count_observables`), even if they cancel out below.
            for k in detectors
                st.maxdet = max(st.maxdet, k)
            end
            for k in instruction.logicals
                st.maxobs = max(st.maxobs, k)
            end
            push!(st.mechanisms, (instruction.p, _odd_parity(detectors), _odd_parity(instruction.logicals)))
        elseif instruction isa DemRepeatInstruction
            for _ in 1:instruction.count
                _dem_interpret!(st, instruction.body)
            end
        end
    end
    return st
end

##############################
# User-facing API
##############################

"""
$(TYPEDSIGNATURES)

Parse the contents of a Stim detector error model (`.dem`) file into a
[`DetectorErrorModelCircuit`](@ref) that can be sampled with [`pftrajectories`](@ref).
See [`read_detector_error_model`](@ref) for reading directly from a file
and for a description of the format and of the returned circuit.

```jldoctest
julia> circuit = parse_detector_error_model(\"\"\"
           detector D0
           detector D1
           logical_observable L0
           error(0.1) D0
           error(0.2) D0 D1 L0
           \"\"\");

julia> circuit.n_detectors, circuit.n_logicals
(2, 1)

julia> length(circuit) # one `DemDeclaration` and two `DetectorError`s
3
```
"""
function parse_detector_error_model(content::AbstractString)
    instructions, _ = _dem_parse_block(split(content, '\n'), 1, false)
    st = DemInterpreterState(0, -1, -1, Tuple{Float64,Vector{Int},Vector{Int}}[])
    _dem_interpret!(st, instructions)
    n_detectors = st.maxdet + 1
    n_logicals = st.maxobs + 1
    ops = AbstractOperation[DemDeclaration(n_detectors, n_logicals)]
    for (p, detectors, logicals) in st.mechanisms
        push!(ops, DetectorError(p, detectors .+ 1, logicals .+ (n_detectors + 1)))
    end
    return DetectorErrorModelCircuit(ops, n_detectors, n_logicals)
end

"""
$(TYPEDSIGNATURES)

Read a Stim detector error model (`.dem`) file and return a
[`DetectorErrorModelCircuit`](@ref): a vector of operations (one no-op
[`DemDeclaration`](@ref) followed by one [`DetectorError`](@ref) per error mechanism)
that is accepted directly by the existing [`pftrajectories`](@ref) methods.

```julia
circuit = read_detector_error_model("surface_code.dem")
frames = pftrajectories(circuit; trajectories=10_000)
detectors = detectorview(circuit, frames)   # trajectories × n_detectors Bool array
logicals = observableview(circuit, frames)  # trajectories × n_logicals Bool array
```

`measurements(frames)` is a `trajectories × (n_detectors + n_logicals)` `Bool` matrix:
detector `Dk` is stored in column `k+1` and logical observable `Lk` in column
`n_detectors + k + 1`; [`detectorview`](@ref) and [`observableview`](@ref) slice the
two blocks apart. Each `error(p)` instruction of the file becomes one independent
[`DetectorError`](@ref) Bernoulli mechanism, so sampling matches the semantics of
Stim's `.dem` documentation (`repeat` blocks are unrolled, `shift_detectors` offsets
are resolved at import time, and targets repeated within one instruction cancel out).

Supported syntax: `error(p) D... L...` (including `^` separators), `detector`,
`logical_observable`, `shift_detectors`, nested `repeat` blocks, instruction tags,
comments, and blank lines. Coordinate arguments and instruction tags are parsed and
discarded; the suggested decompositions marked by `^` separators do not affect
sampling. Unsupported or malformed syntax raises an `ArgumentError` naming
the offending line.

A method accepting an `IO` stream is also provided.

See also: [`parse_detector_error_model`](@ref), [`pftrajectories`](@ref), [`measurements`](@ref)
"""
read_detector_error_model(path::AbstractString) = open(io -> parse_detector_error_model(read(io, String)), path)
read_detector_error_model(io::IO) = parse_detector_error_model(read(io, String))
