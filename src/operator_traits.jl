abstract type OperatorDeterminismTrait end

struct DeterministicOperatorTrait <: OperatorDeterminismTrait end
struct NondeterministicOperatorTrait <: OperatorDeterminismTrait end

"""$(TYPEDSIGNATURES)"""
operatordeterminism(::Type{<:AbstractCliffordOperator}) = DeterministicOperatorTrait()

"""$(TYPEDSIGNATURES)"""
operatordeterminism(::Type{<:AbstractOperation}) = NondeterministicOperatorTrait()
