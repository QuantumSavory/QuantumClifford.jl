"""A single independent error mechanism in a Stim detector error model.

The `detectors` and `logical_observables` fields use Stim's zero-based target
numbering. Thus QuantumClifford stabilizer row `i` is written as detector
target `D(i - 1)`, and fault-matrix row `j` is written as logical target
`L(j - 1)`.
"""
struct DetectorErrorModelError
    probability::Float64
    detectors::Vector{Int}
    logical_observables::Vector{Int}
end

"""A code-capacity Stim detector error model exported from an ECC code.

`num_detectors` is the number of stabilizer checks declared as `D` targets.
`num_observables` is the number of logical-observable rows from
[`faults_matrix`](@ref), declared as `L` targets.
"""
struct DetectorErrorModel
    num_detectors::Int
    num_observables::Int
    errors::Vector{DetectorErrorModelError}
end

function _dem_probability(p, name::Symbol)
    p isa Real || throw(ArgumentError("`$name` must be a real probability."))
    prob = Float64(p)
    isfinite(prob) && 0 <= prob <= 1 ||
        throw(ArgumentError("`$name` must be between 0 and 1."))
    return prob
end

_zero_based_true_indices(bits) = findall(bit -> bit != 0, bits) .- 1

function _y_logical_targets(faults::AbstractMatrix{Bool}, n::Int, qubit::Int)
    targets = Int[]
    for row in axes(faults, 1)
        xor(faults[row, qubit], faults[row, n + qubit]) && push!(targets, row - 1)
    end
    return targets
end

function _push_dem_error!(errors, probability, detectors, logical_observables)
    probability == 0 && return errors
    push!(
        errors,
        DetectorErrorModelError(probability, collect(detectors), collect(logical_observables)),
    )
    return errors
end

"""Build a code-capacity Stim detector error model for an ECC code.

Each stabilizer row in `parity_checks(code)` is declared as a detector target.
Detector targets are zero-based and follow QuantumClifford's stabilizer row
order: stabilizer row `i` maps to detector target `D(i - 1)`.

Logical observable targets are also zero-based and follow the row order of
[`faults_matrix`](@ref): the first `k` rows are logical-X observables, and the
next `k` rows are logical-Z observables. A physical `Y` fault flips the xor of
the corresponding physical `X` and `Z` fault-matrix columns.

Only nonzero independent physical fault probabilities are emitted. Error
mechanisms are ordered deterministically by physical qubit, then by Pauli fault
`X`, `Y`, `Z`.
"""
function detector_error_model(code; px=0.0, py=0.0, pz=0.0)
    px = _dem_probability(px, :px)
    py = _dem_probability(py, :py)
    pz = _dem_probability(pz, :pz)

    checks = parity_checks(code)
    faults = faults_matrix(checks)
    n = code_n(checks)
    errors = DetectorErrorModelError[]

    for qubit in 1:n
        if px != 0
            _push_dem_error!(
                errors,
                px,
                _zero_based_true_indices(comm(checks, single_x(n, qubit))),
                _zero_based_true_indices(@view faults[:, qubit]),
            )
        end
        if py != 0
            _push_dem_error!(
                errors,
                py,
                _zero_based_true_indices(comm(checks, single_y(n, qubit))),
                _y_logical_targets(faults, n, qubit),
            )
        end
        if pz != 0
            _push_dem_error!(
                errors,
                pz,
                _zero_based_true_indices(comm(checks, single_z(n, qubit))),
                _zero_based_true_indices(@view faults[:, n + qubit]),
            )
        end
    end

    return DetectorErrorModel(code_s(checks), size(faults, 1), errors)
end

function _write_dem_targets(io::IO, prefix::String, targets)
    for target in targets
        print(io, " ", prefix, target)
    end
end

"""Write a [`DetectorErrorModel`](@ref) using Stim's `.dem` text format."""
function write_detector_error_model(io::IO, dem::DetectorErrorModel)
    for detector in 0:dem.num_detectors - 1
        println(io, "detector D", detector)
    end
    for observable in 0:dem.num_observables - 1
        println(io, "logical_observable L", observable)
    end
    for err in dem.errors
        print(io, "error(", err.probability, ")")
        _write_dem_targets(io, "D", err.detectors)
        _write_dem_targets(io, "L", err.logical_observables)
        println(io)
    end
    return nothing
end

"""Write a [`DetectorErrorModel`](@ref) to a `.dem` file."""
function write_detector_error_model(path::AbstractString, dem::DetectorErrorModel)
    open(path, "w") do io
        write_detector_error_model(io, dem)
    end
    return nothing
end
