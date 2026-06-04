struct DetectorErrorModelTerm
    probability::Float64
    detectors::Vector{Int}
    observables::Vector{Int}
end

"""A Stim-compatible detector error model for an ECC code.

Detectors are numbered from zero following the row order of `parity_checks(code)`.
Logical observables are numbered from zero following the row order of
`faults_matrix(code)`, i.e. logical-X observables followed by logical-Z
observables.

Use [`detector_error_model`](@ref) to construct a code-capacity model and
[`write_detector_error_model`](@ref) to write Stim `.dem` text.
"""
struct DetectorErrorModel
    num_detectors::Int
    num_observables::Int
    errors::Vector{DetectorErrorModelTerm}
end

Base.length(dem::DetectorErrorModel) = length(dem.errors)
Base.isempty(dem::DetectorErrorModel) = isempty(dem.errors)

function _validate_dem_probability(p::Real, name::Symbol)
    0 <= p <= 1 || throw(DomainError(p, "`$name` must be in the range [0, 1]"))
    return Float64(p)
end

_target_indices(bits) = [i - 1 for i in findall(!iszero, bits)]

function _logical_targets(fm, q::Int, n::Int, pauli::Symbol)
    if pauli === :X
        return _target_indices(@view fm[:, q])
    elseif pauli === :Z
        return _target_indices(@view fm[:, n + q])
    else
        return _target_indices(xor.(@view(fm[:, q]), @view(fm[:, n + q])))
    end
end

"""Create a code-capacity Stim detector error model for `code`.

For each physical qubit, the model includes independent single-qubit `X`, `Y`,
and `Z` Pauli error mechanisms with probabilities `px`, `py`, and `pz`.
Detector targets come from the checks flipped by `comm(error, parity_checks(code))`.
Logical-observable targets come from [`faults_matrix`](@ref), whose rows are
ordered as logical-X observables followed by logical-Z observables.
"""
function detector_error_model(code; px::Real=0, py::Real=0, pz::Real=0)
    px = _validate_dem_probability(px, :px)
    py = _validate_dem_probability(py, :py)
    pz = _validate_dem_probability(pz, :pz)

    H = parity_checks(code)
    s, n = size(H)
    fm = faults_matrix(H)
    size(fm, 2) == 2n || throw(DimensionMismatch("expected faults_matrix with $(2n) columns, got $(size(fm, 2))"))

    terms = DetectorErrorModelTerm[]
    for q in 1:n
        for (probability, pauli, error) in (
            (px, :X, single_x(n, q)),
            (py, :Y, single_y(n, q)),
            (pz, :Z, single_z(n, q)),
        )
            iszero(probability) && continue
            detectors = _target_indices(comm(error, H))
            observables = _logical_targets(fm, q, n, pauli)
            (isempty(detectors) && isempty(observables)) && continue
            push!(terms, DetectorErrorModelTerm(probability, detectors, observables))
        end
    end

    return DetectorErrorModel(s, size(fm, 1), terms)
end

function _print_dem_targets(io::IO, detectors, observables)
    for detector in detectors
        print(io, " D", detector)
    end
    for observable in observables
        print(io, " L", observable)
    end
end

"""Write `dem` as Stim detector error model text to `io`.

The writer emits deterministic `detector` and `logical_observable` declarations
before the `error(p) ...` lines.
"""
function write_detector_error_model(io::IO, dem::DetectorErrorModel)
    for detector in 1:dem.num_detectors
        println(io, "detector D", detector - 1)
    end
    for observable in 1:dem.num_observables
        println(io, "logical_observable L", observable - 1)
    end
    for term in dem.errors
        print(io, "error(", repr(term.probability), ")")
        _print_dem_targets(io, term.detectors, term.observables)
        println(io)
    end
    return nothing
end

function write_detector_error_model(path::AbstractString, dem::DetectorErrorModel)
    open(path, "w") do io
        write_detector_error_model(io, dem)
    end
    return nothing
end

function write_detector_error_model(io::IO, code; kwargs...)
    return write_detector_error_model(io, detector_error_model(code; kwargs...))
end

function write_detector_error_model(path::AbstractString, code; kwargs...)
    return write_detector_error_model(path, detector_error_model(code; kwargs...))
end
