using Random: AbstractRNG, GLOBAL_RNG

abstract type AbstractSymbolicOperator <: AbstractCliffordOperator end
abstract type AbstractSingleQubitOperator <: AbstractSymbolicOperator end
abstract type AbstractTwoQubitOperator <: AbstractSymbolicOperator end

# Stim has a good list of specialized single and two qubit operations at https://github.com/quantumlib/Stim/blob/e51ea66d213b25920e72c08e53266ec56fd14db4/src/stim/stabilizers/tableau_specialized_prepend.cc
# Note that their specialized operations are for prepends (right multiplications), while we implement append (left multiplication) operations.

"""A "symbolic" Hadamard operator which permits faster multiplication than an operator expressed as an explicit tableau.

```jldoctest
julia> sh = sHadamard(1)
Symbolic single-qubit gate on qubit 1
X ⟼ + Z
Z ⟼ + X

julia> typeof(sh)
sHadamard

julia> h = CliffordOperator(sh,1) # Transforming it back into an explicit tableau representation
X ⟼ + Z
Z ⟼ + X

julia> typeof(h)
CliffordOperator{Vector{UInt8}, Matrix{UInt64}}
```

See also: [`SingleQubitOperator`](@ref)
"""
struct sHadamard <: AbstractSingleQubitOperator
    q::Int
end

"""A "symbolic" Phase operator which permits faster multiplication than an operator expressed as an explicit tableau.

```jldoctest
julia> sp = sPhase(1)
Symbolic single-qubit gate on qubit 1
X ⟼ + Y
Z ⟼ + Z

julia> typeof(sp)
sPhase

julia> p = CliffordOperator(sp,1) # Transforming it back into an explicit tableau representation
X ⟼ + Y
Z ⟼ + Z

julia> typeof(p)
CliffordOperator{Vector{UInt8}, Matrix{UInt64}}
```

See also: [`SingleQubitOperator`](@ref)
"""
struct sPhase <: AbstractSingleQubitOperator
    q::Int
end

"""A "symbolic" inverse of the Phase operator.

See also: [`sPhase`](@ref), [`SingleQubitOperator`](@ref)
"""
struct sInvPhase <: AbstractSingleQubitOperator
    q::Int
end

"""A "symbolic" single-qubit Identity operation.

See also: [`SingleQubitOperator`](@ref)
"""
struct sId1 <: AbstractSingleQubitOperator
    q::Int
end

"""A "symbolic" single-qubit Pauli X operation.

See also: [`SingleQubitOperator`](@ref)
"""
struct sX <: AbstractSingleQubitOperator
    q::Int
end
"""A "symbolic" single-qubit Pauli Y operation.

See also: [`SingleQubitOperator`](@ref)
"""
struct sY <: AbstractSingleQubitOperator
    q::Int
end
"""A "symbolic" single-qubit Pauli Z operation.

See also: [`SingleQubitOperator`](@ref)
"""
struct sZ <: AbstractSingleQubitOperator
    q::Int
end

"""A "symbolic" general single-qubit operator which permits faster multiplication than an operator expressed as an explicit tableau.

```jldoctest
julia> op = SingleQubitOperator(2, true, true, true, false, true, true) # Tableau components and phases
Symbolic single-qubit gate on qubit 2
X ⟼ - Y
Z ⟼ - X

julia> typeof(op)
SingleQubitOperator

julia> t_op = CliffordOperator(op, 3) # Transforming it back into an explicit tableau representation (specifying the size)
X__ ⟼ + X__
_X_ ⟼ - _Y_
__X ⟼ + __X
Z__ ⟼ + Z__
_Z_ ⟼ - _X_
__Z ⟼ + __Z

julia> typeof(t_op)
CliffordOperator{Vector{UInt8}, Matrix{UInt64}}

julia> CliffordOperator(op, 1, compact=true) # You can also extract just the non-trivial part of the tableau
X ⟼ - Y
Z ⟼ - X
```

See also: [`sHadamard`](@ref), [`sPhase`](@ref), [`sId1`](@ref), [`sX`](@ref), [`sY`](@ref), [`sZ`](@ref), [`CliffordOperator`](@ref)

Or simply consult `subtypes(QuantumClifford.AbstractSingleQubitOperator)` and
`subtypes(QuantumClifford.AbstractTwoQubitOperator)` for a full list.
You can think of the `s` prefix as "symbolic" or "sparse".
"""
struct SingleQubitOperator <: AbstractSingleQubitOperator
    q::Int
    xx::Bool
    xz::Bool
    zx::Bool
    zz::Bool
    px::Bool
    pz::Bool
end

SingleQubitOperator(h::sHadamard) = SingleQubitOperator(h.q, false, true , true , false, false, false)
SingleQubitOperator(p::sPhase)    = SingleQubitOperator(p.q, true , true , false, true , false, false)
SingleQubitOperator(p::sInvPhase) = SingleQubitOperator(p.q, true , true , false, true , true , false)
SingleQubitOperator(p::sId1)      = SingleQubitOperator(p.q, true , false, false, true , false, false)
SingleQubitOperator(p::sX)        = SingleQubitOperator(p.q, true , false, false, true , false, true)
SingleQubitOperator(p::sY)        = SingleQubitOperator(p.q, true , false, false, true , true , true)
SingleQubitOperator(p::sZ)        = SingleQubitOperator(p.q, true , false, false, true , true , false)
SingleQubitOperator(o::SingleQubitOperator) = o
function SingleQubitOperator(op::CliffordOperator, qubit)
    nqubits(op)==1 || throw(DimensionMismatch("You are trying to convert a multiqubit `CliffordOperator` into a symbolic `SingleQubitOperator`."))
    SingleQubitOperator(qubit,op.tab[1,1]...,op.tab[2,1]...,(~).(iszero.(op.tab.phases))...)
end
SingleQubitOperator(op::CliffordOperator) = SingleQubitOperator(op, 1)
CliffordOperator(op::AbstractSingleQubitOperator, n; kw...) = CliffordOperator(SingleQubitOperator(op), n; kw...)
function CliffordOperator(op::SingleQubitOperator, n; compact=false)
    if compact
        n==1 || throw(ArgumentError("Set `n=1` as a `SingleQubitOperator` being compacted (`compact=true`) has to result in a 1D `CliffordOperator`."))
        return CliffordOperator(Stabilizer([op.px ? 0x2 : 0x0, op.pz ? 0x2 : 0x0],[op.xx op.xz; op.zx op.zz]))
    else
        n >= op.q || throw(DimensionMismatch("Set a larger `n`, otherwise the `SingleQubitOperator` can not fit in the allocated `CliffordOperator`. Use `compact=true` if you want to discard the target index."))
        c = one(CliffordOperator, n)
        c[op.q,op.q] = op.xx, op.xz # TODO define an `embed` helper function
        c[n+op.q,op.q] = op.zx, op.zz
        c.tab.phases[op.q] = op.px ? 0x2 : 0x0 # TODO define a `phasesview` or `phases` helper function 
        c.tab.phases[n+op.q] = op.pz ? 0x2 : 0x0
        return c
    end
end
function Base.show(io::IO, op::AbstractSingleQubitOperator)
    print(io, "Symbolic single-qubit gate on qubit $(op.q)\n")
    show(io, CliffordOperator(op, 1, compact=true))
end

#TODO these are too low level to be used directly in the apply! functions. Make a better abstractions.
@inline getshift(Tme::Type,col::Int) = _mod(Tme,col-1)
@inline getmask(Tme::Type,col::Int) = Tme(0x1)<<getshift(Tme,col)
@inline getbigindex(Tme::Type,col::Int) = _div(Tme,col-1)+1

Base.@propagate_inbounds function getxbit(s, r, c)
    Tme = eltype(s.xzs)
    s.xzs[getbigindex(Tme,c),r]&getmask(Tme,c)
end
Base.@propagate_inbounds function getzbit(s, r, c)
    Tme = eltype(s.xzs)
    s.xzs[end÷2+getbigindex(Tme,c),r]&getmask(Tme,c)
end
Base.@propagate_inbounds function setxbit(s, r, c, x)
    Tme = eltype(s.xzs)
    cbig = getbigindex(Tme,c)
    s.xzs[cbig,r] &= ~getmask(Tme,c)
    s.xzs[cbig,r] |= x
end
Base.@propagate_inbounds function setzbit(s, r, c, z)
    Tme = eltype(s.xzs)
    cbig = getbigindex(Tme,c)
    s.xzs[end÷2+cbig,r] &= ~getmask(Tme,c)
    s.xzs[end÷2+cbig,r] |= z
end
Base.@propagate_inbounds setxbit(s, r, c, x, shift) = setxbit(s, r, c, x<<shift)
Base.@propagate_inbounds setzbit(s, r, c, z, shift) = setzbit(s, r, c, z<<shift)

function _apply!(stab::AbstractStabilizer, op::AbstractSymbolicOperator; phases::Bool=true)
    apply!(stab,op; phases=Val(phases))
end

function _apply!(stab::AbstractStabilizer, h::sHadamard; phases::Val{B}=Val(true)) where B
    s = tab(stab)
    c = h.q
    @batch per=core minbatch=200 for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        setxbit(s, r, c, z)
        setzbit(s, r, c, x)
        B && x!=0 && z!=0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function _apply!(stab::AbstractStabilizer, p::sPhase; phases::Val{B}=Val(true)) where B
    s = tab(stab)
    c = p.q
    @batch per=core minbatch=200 for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        #setxbit no op
        setzbit(s, r, c, x⊻z)
        B && x!=0 && z!=0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function _apply!(stab::AbstractStabilizer, p::sInvPhase; phases::Val{B}=Val(true)) where B
    s = tab(stab)
    c = p.q
    @batch per=core minbatch=200 for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        #setxbit no op
        setzbit(s, r, c, x⊻z)
        B && x!=0 && z==0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function _apply!(stab::AbstractStabilizer, p::sX; phases::Val{B}=Val(true)) where B
    B || return stab
    s = tab(stab)
    c = p.q
    @batch per=core minbatch=200 for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        B && z!=0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function _apply!(stab::AbstractStabilizer, p::sZ; phases::Val{B}=Val(true)) where B
    B || return stab
    s = tab(stab)
    c = p.q
    @batch per=core minbatch=200 for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        B && x!=0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function _apply!(stab::AbstractStabilizer, p::sY; phases::Val{B}=Val(true)) where B
    B || return stab
    s = tab(stab)
    c = p.q
    @batch per=core minbatch=200 for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        B && (x⊻z)!=0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function _apply!(stab::AbstractStabilizer, ::sId1; phases::Val{B}=Val(true)) where B
    stab
end

function _apply!(stab::AbstractStabilizer, op::SingleQubitOperator; phases::Val{B}=Val(true)) where B # TODO Generated functions that simplify the whole `if phases` branch might be a good optimization, but a quick benchmakr comparing sHadamard to SingleQubitOperator(sHadamard) did not show a worthwhile difference.
    s = tab(stab)
    c = op.q
    Tme = eltype(s.xzs)
    sh = getshift(Tme, c)
    xx,zx,xz,zz = Tme.((op.xx,op.zx,op.xz,op.zz)) .<< sh
    anticom = ~iszero((~zz & xz & ~xx & zx) | ( zz & ~xz & xx & zx) | (zz &  xz & xx & ~zx))
    @batch per=core minbatch=200 for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        setxbit(s, r, c, (x&xx)⊻(z&zx))
        setzbit(s, r, c, (x&xz)⊻(z&zz))
        if B
            if op.px && ~iszero(x)
                s.phases[r] += 0x2
                s.phases[r] &= 3
            end
            if op.pz && ~iszero(z)
                s.phases[r] += 0x2
                s.phases[r] &= 3
            end
            if ~iszero(x&z) && anticom
                s.phases[r] += 0x2
                s.phases[r] &= 3
            end
        end
    end
    stab
end

const all_single_qubit_patterns = (
    (true, false, false, true), # X, Z ↦ X, Z
    (false, true, true, true),  # X, Z ↦ Z, Y
    (true, true, true, false),  # X, Z ↦ Y, X
    (false, true, true, false), # X, Z ↦ Z, X - Hadamard
    (true, false, true, true),  # X, Z ↦ X, Y
    (true, true, false, true)   # X, Z ↦ Y, Z - Phase
)

"""Generate a symbolic single-qubit gate given its index. Optionally, set non-trivial phases."""
function enumerate_single_qubit_gates(index; qubit=1, phases=(false,false))
    @assert index<=6 "Only 6 single-qubit gates exit, up to the choice of phases"
    if phases==(false,false)
        if index==4
            return sHadamard(qubit)
        elseif index==6
            return sPhase(qubit)
        else
            return SingleQubitOperator(qubit, all_single_qubit_patterns[index]..., false, false)
        end
    else
        if index==1
            if     (phases[1], phases[2]) == (false, false)
                return sId1(qubit)
            elseif (phases[1], phases[2]) == (false,  true)
                return sX(qubit)
            elseif (phases[1], phases[2]) == (true,  false)
                return sZ(qubit)
            else
                return sY(qubit)
            end
        else
            return SingleQubitOperator(qubit, all_single_qubit_patterns[index]..., phases...)
        end
    end
end

"""Random symbolic single-qubit Clifford applied to qubit at index `qubit`.

See also: [`SingleQubitOperator`](@ref), [`random_clifford`](@ref)
"""
function random_clifford1(rng::AbstractRNG, qubit)
    return enumerate_single_qubit_gates(rand(rng,1:6),qubit=qubit,phases=rand(rng,Bool,2))
end
random_clifford1(qubit) = random_clifford1(GLOBAL_RNG, qubit)


"""A "symbolic" CNOT

See also: [`SingleQubitOperator`](@ref)
"""
struct sCNOT <: AbstractTwoQubitOperator
    q1::Int
    q2::Int
end

"""A "symbolic" SWAP

See also: [`SingleQubitOperator`](@ref)
"""
struct sSWAP <: AbstractTwoQubitOperator
    q1::Int
    q2::Int
end

function _apply!(stab::AbstractStabilizer, cnot::sCNOT; phases::Val{B}=Val(true)) where B
    s = tab(stab)
    q1 = cnot.q1
    q2 = cnot.q2
    Tme = eltype(s.xzs)
    shift = getshift(Tme, q1) - getshift(Tme, q2)
    @batch per=core minbatch=200 for r in eachindex(s)
        x1 = getxbit(s, r, q1)
        z1 = getzbit(s, r, q1)
        x2 = getxbit(s, r, q2)
        z2 = getzbit(s, r, q2)
        setxbit(s, r, q2, x2⊻(x1<<-shift))
        setzbit(s, r, q1, z1⊻(z2<<shift))
        if B && ~iszero( x1&z1&((x2&z2)<<shift) | (x1&z2<<shift & ~(z1|x2<<shift)) )
            s.phases[r] += 0x2
            s.phases[r] &= 3
        end
    end
    stab
end

function _apply!(stab::AbstractStabilizer, swap::sSWAP; phases::Val{B}=Val(true)) where B
    s = tab(stab)
    q1 = swap.q1
    q2 = swap.q2
    Tme = eltype(s.xzs)
    shift = getshift(Tme, q1) - getshift(Tme, q2)
    @batch per=core minbatch=200 for r in eachindex(s)
        x1 = getxbit(s, r, q1)
        z1 = getzbit(s, r, q1)
        x2 = getxbit(s, r, q2)
        z2 = getzbit(s, r, q2)
        setxbit(s, r, q1, x2, shift)
        setzbit(s, r, q1, z2, shift)
        setxbit(s, r, q2, x1, -shift)
        setzbit(s, r, q2, z1, -shift)
    end
    stab
end

##############################
# Functions that perform direct application of common operators without needing an operator instance
##############################
# TODO is there a reason to keep these given that sZ/sX/sY exist?
# TODO currently these are faster than sZ/sX/sY

"""Apply a Pauli Z to the `i`-th qubit of state `s`. You should use `apply!(stab,sZ(i))` instead of this."""
function apply_single_z!(stab::AbstractStabilizer, i)
    s = tab(stab)
    Tme = eltype(s.xzs)
    bigi = _div(Tme,i-1)+1
    smalli = _mod(Tme,i-1)
    mask = Tme(0x1)<<smalli
    @inbounds @simd for row in 1:size(s.xzs,2)
        if !iszero(s.xzs[bigi,row] & mask)
            s.phases[row] = (s.phases[row]+0x2)&0x3
        end
    end
    stab
end

"""Apply a Pauli X to the `i`-th qubit of state `s`. You should use `apply!(stab,sX(i))` instead of this."""
function apply_single_x!(stab::AbstractStabilizer, i)
    s = tab(stab)
    Tme = eltype(s.xzs)
    bigi = _div(Tme,i-1)+1
    smalli = _mod(Tme,i-1)
    mask = Tme(0x1)<<smalli
    @inbounds @simd for row in 1:size(s.xzs,2)
        if !iszero(s.xzs[end÷2+bigi,row] & mask)
            s.phases[row] = (s.phases[row]+0x2)&0x3
        end
    end
    stab
end

"""Apply a Pauli Y to the `i`-th qubit of state `s`. You should use `apply!(stab,sY(i))` instead of this."""
function apply_single_y!(stab::AbstractStabilizer, i)
    s = tab(stab)
    Tme = eltype(s.xzs)
    bigi = _div(Tme,i-1)+1
    smalli = _mod(Tme,i-1)
    mask = Tme(0x1)<<smalli
    @inbounds @simd for row in 1:size(s.xzs,2)
        if !iszero((s.xzs[bigi,row] & mask) ⊻ (s.xzs[end÷2+bigi,row] & mask))
            s.phases[row] = (s.phases[row]+0x2)&0x3
        end
    end
    stab
end
