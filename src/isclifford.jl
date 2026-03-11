"""
$(SIGNATURES)

Trait function to determine if an operation is a Clifford gate.
Users can extend this for custom gate types by defining new methods.

# Examples
```jldoctest
julia> isclifford(sHadamard(1))
true

julia> isclifford(sPhase(1))
true

julia> isclifford(sCNOT(1,2))
true
```
"""
function isclifford end

isclifford(::AbstractCliffordOperator) = true
isclifford(::AbstractSingleQubitOperator) = true
isclifford(::AbstractTwoQubitOperator) = true
isclifford(::AbstractSymbolicOperator) = true

isclifford(::sMX) = true
isclifford(::sMY) = true
isclifford(::sMZ) = true
isclifford(::sMRX) = true
isclifford(::sMRY) = true
isclifford(::sMRZ) = true
isclifford(::PauliMeasurement) = true

isclifford(op::SparseGate) = isclifford(op.cliff)

isclifford(::AbstractOperation) = false
