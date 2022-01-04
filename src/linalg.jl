"""The inner product of two Stabilizers.

Based on [garcia2012efficient](@cite).

See also: [`logdot`](@ref)"""
function LinearAlgebra.dot(s1::AbstractStabilizer, s2::AbstractStabilizer)
    ld = logdot(s1,s2)
    if isnothing(ld)
        return 0.0
    else
        return 2.0^(-ld/2)
    end
end

"""Logarithm of the inner product between to Stabilizer states.

If the result is `nothing`, the dot inner product is zero.
Otherwise the inner product is `2^(-logdot/2)`.

The actual inner product can be computed with `LinearAlgebra.dot`.

Based on [garcia2012efficient](@cite)."""
function logdot(s1::AbstractStabilizer, s2::AbstractStabilizer) # TODO verify rank
    logdot(stabilizerview(s1),stabilizerview(s2))
end

function logdot(s1::Stabilizer, s2::Stabilizer)
    if nqubits(s1)!=length(s1) || nqubits(s2)!=length(s2) # TODO implement this
        throw(DomainError("Only pure (not mixed) states are supported when calculating inner product."))
    end
    if nqubits(s1)!=nqubits(s2)
        throw(DimensionMismatch("Inner product can be calculated only between states with the same number of qubits."))
    end
    c1_inv = inv(CliffordOperator(copy(s1)))
    s2_prime = canonicalize!(c1_inv*s2)
    canonicalize!(s2_prime)
    k = 0
    for row in eachindex(s2_prime)
        if any(i->s2_prime[row,i][1], 1:nqubits(s2_prime)) # X or Y
            k += 1
        else
            if !iszero(s2_prime.phases[row])
                return nothing
            end
        end
    end
    return k
end

LinearAlgebra.rank(s::Stabilizer)   = throw(BadDataStructure("Using a `Stabilizer` type does not permit automatic tracking of the rank. Use `MixedDestabilizer` type instead or track the rank manually.",
                                            :rank, :Stabilizer))

LinearAlgebra.rank(s::Destabilizer) = throw(BadDataStructure("Using a `Destabilizer` type does not permit automatic tracking of the rank. Use `MixedDestabilizer` type instead or track the rank manually.",
                                            :rank, :Stabilizer))

LinearAlgebra.rank(s::MixedStabilizer) = s.rank
LinearAlgebra.rank(s::MixedDestabilizer) = s.rank