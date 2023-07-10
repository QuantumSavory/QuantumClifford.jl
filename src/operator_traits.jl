abstract type OperatorDeterminismTrait end

struct DeterministicOperatorTrait <: OperatorDeterminismTrait end
struct NondeterministicOperatorTrait <: OperatorDeterminismTrait end

operatordeterminism(::Type{<:AbstractCliffordOperator}) = DeterministicOperatorTrait()
operatordeterminism(::Type{<:AbstractOperation}) = NondeterministicOperatorTrait()
