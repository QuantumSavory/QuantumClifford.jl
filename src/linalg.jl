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

LinearAlgebra.rank(s::Stabilizer)   = throw(BadDataStructure("Using a `Stabilizer` type does not permit automatic tracking of the rank. Use `length`, the `MixedDestabilizer` type, or track the rank manually.",
                                            :rank, :Stabilizer))

LinearAlgebra.rank(s::Destabilizer) = throw(BadDataStructure("Using a `Destabilizer` type does not permit automatic tracking of the rank. Use `length`, the `MixedDestabilizer` type, or track the rank manually.",
                                            :rank, :Stabilizer))

LinearAlgebra.rank(s::MixedStabilizer) = s.rank
LinearAlgebra.rank(s::MixedDestabilizer) = s.rank

"""Tensor product between operators or tableaux. See also [`tensor`](@ref) and [`tensor_pow`](@ref)."""
function ⊗ end

function ⊗(ops::AbstractStabilizer...) # TODO optimize this by doing conversion to common type to enable preallocation
    foldl(⊗, ops[2:end], init=ops[1])
end

"""Tensor product between operators or tableaux. See also [`⊗`](@ref) and [`tensor_pow`](@ref)."""
const tensor = ⊗

"""Repeated tensor product of an operators or a tableau. See also [`⊗`](@ref) and [`tensor_pow`](@ref)."""
function tensor_pow(op,power)
    if power==1
        return op
    else
        return tensor((op for i in 1:power)...)
    end
end

function ⊗(ops::Stabilizer...)
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

function ⊗(ops::MixedDestabilizer...)
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops)
    rtot = sum(LinearAlgebra.rank, ops)
    tab = zero(Stabilizer, 2*ntot, ntot)
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

⊗(ops::AbstractCliffordOperator...) = ⊗(CliffordOperator.(ops)...)

function ⊗(ops::CliffordOperator...) # TODO implement \otimes for Destabilizer and use it here
    length(ops)==1 && return ops[1]
    ntot = sum(nqubits, ops)
    tab = zero(Stabilizer, 2*ntot, ntot)
    last_zrow = ntot
    last_xrow = 0
    for op in ops
        t = op.tab
        _, last_zrow, _ = puttableau!(tab, (@view t[end÷2+1:end]), last_zrow, last_xrow)
        _, last_xrow, _ = puttableau!(tab, (@view t[1:end÷2]), last_xrow, last_xrow)
    end
    CliffordOperator(tab)
end

"""Put source tableau in target tableau at given row and column. Assumes target location is zeroed out."""
@inline function puttableau!(target::Stabilizer{V1,M1}, source::Stabilizer{V2,M2}, row::Int, col::Int; phases::Bool=true) where {V1,V2,T<:Unsigned,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T}}
    xzs = target.xzs
    ph = target.phases
    sxzs = source.xzs
    sph = source.phases
    r,n = size(source)
    bₗ = _div(T,col)+1
    bᵣ = bₗ + 1
    eₗ = bₗ + _div(T,n-1)
    eᵣ = _div(T,col+n-1)+1
    shiftₗ = _mod(T,col)
    shiftᵣ = 8*sizeof(T)-shiftₗ        
    for i in 1:r
    @inbounds @simd for j in 0:eₗ-bₗ
        xzs[bₗ+j,row+i] |= sxzs[j+1,i] >>> -shiftₗ
        xzs[end÷2+bₗ+j,row+i] |= sxzs[end÷2+j+1,i] >>> -shiftₗ
    end
    @inbounds @simd for j in 0:eᵣ-bᵣ
        xzs[bᵣ+j,row+i] |= sxzs[j+1,i] >>> shiftᵣ
        xzs[end÷2+bᵣ+j,row+i] |= sxzs[end÷2+j+1,i] >>> shiftᵣ
    end
    end
    phases && (ph[row+1:row+r] .= sph)
    target, row+r, col+n
end