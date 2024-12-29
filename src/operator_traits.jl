abstract type OperatorDeterminismTrait end

struct DeterministicOperatorTrait <: OperatorDeterminismTrait end
struct NondeterministicOperatorTrait <: OperatorDeterminismTrait end

operatordeterminism(::Type{<:AbstractCliffordOperator}) = DeterministicOperatorTrait()
operatordeterminism(::Type{<:AbstractOperation}) = NondeterministicOperatorTrait()
operatordeterminism(::Type{sMZ}) = NondeterministicOperatorTrait()
operatordeterminism(::Type{sMX}) = NondeterministicOperatorTrait()
operatordeterminism(::Type{sMY}) = NondeterministicOperatorTrait()