_zero_based_true_indices(bits) = findall(!iszero, bits) .- 1

function _y_logical_targets(faults::AbstractMatrix{Bool}, n::Int, qubit::Int)
    return [
        row - 1 for row in axes(faults, 1) if
        xor(faults[row, qubit], faults[row, n + qubit])
    ]
end

function _write_dem_error(io::IO, probability, detectors, logical_observables)
    probability == 0 && return nothing
    print(io, "error(", probability, ")")
    for detector in detectors
        print(io, " D", detector)
    end
    for observable in logical_observables
        print(io, " L", observable)
    end
    println(io)
    return nothing
end

"""
Build a code-capacity Stim detector error model string for an ECC code.

Each stabilizer row in `parity_checks(code)` is declared as a detector target.
Detector targets are zero-based and follow QuantumClifford's stabilizer row
order: stabilizer row `i` maps to detector target `D(i - 1)`.

Logical observable targets are also zero-based and follow the row order of
[`faults_matrix`](@ref): the first `k` rows are logical-X observables, and the
next `k` rows are logical-Z observables. A physical `Y` fault flips the XOR of
the corresponding physical `X` and `Z` fault-matrix columns.

Only nonzero independent physical fault probabilities are emitted. Error
mechanisms are ordered deterministically by physical qubit, then by Pauli fault
`X`, `Y`, `Z`.

# Example

```julia
using QuantumClifford.ECC

dem = detector_error_model(Steane7(); px=1e-3, pz=1e-3)
write_detector_error_model("steane7.dem", dem)
```
"""
function detector_error_model(code; px=0.0, py=0.0, pz=0.0)
    0 <= px <= 1 && 0 <= py <= 1 && 0 <= pz <= 1 ||
        throw(ArgumentError("`px`, `py`, and `pz` must be between 0 and 1."))
    checks = parity_checks(code)
    faults = faults_matrix(checks)
    n = code_n(checks)
    io = IOBuffer()
    for detector in 0:code_s(checks)-1
        println(io, "detector D", detector)
    end
    for observable in 0:size(faults, 1)-1
        println(io, "logical_observable L", observable)
    end
    for qubit in 1:n
        _write_dem_error(
            io, px,
            _zero_based_true_indices(comm(checks, single_x(n, qubit))),
            _zero_based_true_indices(@view faults[:, qubit]),
        )
        _write_dem_error(
            io, py,
            _zero_based_true_indices(comm(checks, single_y(n, qubit))),
            _y_logical_targets(faults, n, qubit),
        )
        _write_dem_error(
            io, pz,
            _zero_based_true_indices(comm(checks, single_z(n, qubit))),
            _zero_based_true_indices(@view faults[:, n+qubit]),
        )
    end
    return String(take!(io))
end

"""Write a Stim detector error model string to an `IO` object."""
function write_detector_error_model(io::IO, dem::AbstractString)
    write(io, dem)
    return nothing
end

"""Write a Stim detector error model string to a `.dem` file."""
function write_detector_error_model(path::AbstractString, dem::AbstractString)
    open(path, "w") do io
        write_detector_error_model(io, dem)
    end
    return nothing
end
