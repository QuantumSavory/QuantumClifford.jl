
apply!(::Register, ::Any, ::Tuple{}) = error()
apply!(::MixedDestabilizer, ::AbstractCliffordOperator, ::Type{<:AbstractSymbolicOperator}) = error()