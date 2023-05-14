"""
$TYPEDSIGNATURES

Inverse of a `CliffordOperator`
"""
function LinearAlgebra.inv(c::CliffordOperator; phases=true)
    ci = zero(c)
    n = nqubits(c)
    # TODO this transpose can be much faster with proper SIMDing
    for i in 1:n
        for j in 1:n
            ci.tab[i,j] = c.tab[n+j,i][2], c.tab[j,i][2]
            ci.tab[n+i,j] = c.tab[n+j,i][1], c.tab[j,i][1]
        end
    end
    if phases
        ci*c*ci # TODO perform this inplace as in Stim https://github.com/quantumlib/Stim/blob/e51ea66d213b25920e72c08e53266ec56fd14db4/src/stim/stabilizers/tableau.cc#L383
    else
        ci
    end
end

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
function logdot(s1::AbstractStabilizer, s2::AbstractStabilizer) # TODO verify rank # TODO this is currently very inefficient as we discard the destabilizers and then recreate them
    logdot(stabilizerview(s1),stabilizerview(s2))
end

function logdot(s1::Stabilizer, s2::Stabilizer)
    if nqubits(s1)!=length(s1) || nqubits(s2)!=length(s2) # TODO implement this
        throw(DomainError("Only pure (not mixed) states are supported when calculating inner product."))
    end
    if nqubits(s1)!=nqubits(s2)
        throw(DimensionMismatch("Inner product can be calculated only between states with the same number of qubits."))
    end
    c1_inv = inv(CliffordOperator(tab(MixedDestabilizer(copy(s1)))))
    s2_prime = canonicalize!(c1_inv*s2)
    canonicalize!(s2_prime)
    k = 0
    for row in eachindex(s2_prime)
        if any(i->s2_prime[row,i][1], 1:nqubits(s2_prime)) # X or Y
            k += 1
        else
            if !iszero(tab(s2_prime).phases[row])
                return nothing
            end
        end
    end
    return k
end

LinearAlgebra.rank(s::Stabilizer)   = throw(BadDataStructure("Using a `Stabilizer` type does not permit automatic tracking of the rank. Use `length`, `trusted_rank`, the `MixedDestabilizer` type, or track the rank manually.",
                                            :rank, :Stabilizer))
LinearAlgebra.rank(s::Destabilizer) = throw(BadDataStructure("Using a `Destabilizer` type does not permit automatic tracking of the rank. Use `length`, `trusted_rank`, the `MixedDestabilizer` type, or track the rank manually.",
                                            :rank, :Destabilizer))
LinearAlgebra.rank(s::MixedStabilizer) = s.rank
LinearAlgebra.rank(s::MixedDestabilizer) = s.rank

"""A "trusted" `rank` which returns `rank(state)` for `Mixed[De]Stabilizer` and `lenght(state)` for `[De]Stabilizer`."""
function trusted_rank end
trusted_rank(s::Stabilizer) = length(s)
trusted_rank(s::Destabilizer) = length(s)
trusted_rank(s::MixedStabilizer) = LinearAlgebra.rank(s)
trusted_rank(s::MixedDestabilizer) = LinearAlgebra.rank(s)

"""Tensor product between operators or tableaux. See also [`tensor_pow`](@ref)."""
function tensor end

function tensor(ops::AbstractStabilizer...) # TODO optimize this by doing conversion to common type to enable preallocation
    foldl(⊗, ops[2:end], init=ops[1])
end

"""Repeated tensor product of an operators or a tableau. See also [`tensor`](@ref)."""
function tensor_pow(op::Union{<:AbstractStabilizer,<:AbstractCliffordOperator},power)
    if power==1
        return op
    else
        return tensor((op for i in 1:power)...)
    end
end

function tensor(ops::Stabilizer...)
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops)
    rtot = sum(length, ops)
    tab = zero(Stabilizer, rtot, ntot)
    last_row = 0
    last_col = 0
    for op in ops
        _, last_row, last_col = puttableau!(tab, op, last_row, last_col)
    end
    tab
end

function tensor(ops::MixedDestabilizer...)
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops)
    rtot = sum(LinearAlgebra.rank, ops)
    tab = zero(Tableau, 2*ntot, ntot)
    last_svrow = ntot
    last_dvrow = 0
    last_lxrow = rtot
    last_lzrow = ntot+rtot
    last_col = 0
    for op in ops
        _, last_svrow, _        = puttableau!(tab,   stabilizerview(op), last_svrow, last_col)
        _, last_dvrow, _        = puttableau!(tab, destabilizerview(op), last_dvrow, last_col)
        _, last_lxrow, _        = puttableau!(tab,     logicalxview(op), last_lxrow, last_col)
        _, last_lzrow, last_col = puttableau!(tab,     logicalzview(op), last_lzrow, last_col)
    end
    MixedDestabilizer(tab, rtot)
end

tensor(ops::AbstractCliffordOperator...) = ⊗(CliffordOperator.(ops)...)

function tensor(ops::CliffordOperator...) # TODO implement \otimes for Destabilizer and use it here
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops)
    tab = zero(Tableau, 2*ntot, ntot)
    last_zrow = ntot
    last_xrow = 0
    for op in ops
        t = op.tab
        _, last_zrow, _ = puttableau!(tab, (@view t[end÷2+1:end]), last_zrow, last_xrow)
        _, last_xrow, _ = puttableau!(tab, (@view t[1:end÷2]), last_xrow, last_xrow)
    end
    CliffordOperator(tab)
end
