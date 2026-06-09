# TODO: implement apply!(::GeneralizedStabilizer, ::sT)
"""
$(TYPEDEF)

T gate (π/8 phase rotation). This is a non-Clifford gate with optimal
stabilizer extent ξ(T) = (cos(π/8) + tan(π/8)sin(π/8))² ≈ 1.172.

The T gate is diagonal in the computational basis:
T = diag(1, e^(iπ/4)) = R_z(π/4)

$(TYPEDFIELDS)

# Example
```jldoctest
julia> circuit = [sHadamard(1), sT(1), sHadamard(1)];

julia> result = emtrajectories(circuit, 1; trajectories=100);

julia> size(measurements(result))
(100, 1)
```
"""
struct sT <: AbstractNonCliffordOperator
    "Target qubit index (1-indexed)"
    qubit::Int

    function sT(q::Int)
        q > 0 || throw(ArgumentError("Qubit index must be positive, got $q"))
        new(q)
    end
end

nqubits(::sT) = 1

# TODO: implement apply!(::GeneralizedStabilizer, ::sCCZ)
"""
$(TYPEDEF)

Controlled-Controlled-Z gate. This is a non-Clifford gate with optimal
stabilizer extent ξ(CCZ) = 16/9 ≈ 1.778.

The CCZ gate applies a phase of -1 when all three qubits are in state |1⟩:
CCZ|xyz⟩ = (-1)^(xyz)|xyz⟩

$(TYPEDFIELDS)

# Example
```jldoctest
julia> circuit = [sHadamard(1), sHadamard(2), sHadamard(3), sCCZ(1, 2, 3)];

julia> result = emtrajectories(circuit, 3; trajectories=100);

julia> size(measurements(result))
(100, 3)
```
"""
struct sCCZ <: AbstractNonCliffordOperator
    "Target qubit indices (sorted, 1-indexed)"
    qubits::NTuple{3, Int}

    function sCCZ(q1::Int, q2::Int, q3::Int)
        all(q -> q > 0, (q1, q2, q3)) || throw(ArgumentError("All qubit indices must be positive"))
        length(unique((q1, q2, q3))) == 3 || throw(ArgumentError("CCZ requires 3 distinct qubits, got $q1, $q2, $q3"))
        new(Tuple(sort([q1, q2, q3])))
    end
end

sCCZ(qubits::Vector{Int}) = sCCZ(qubits[1], qubits[2], qubits[3])

nqubits(::sCCZ) = 3
