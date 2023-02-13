using Random: AbstractRNG, GLOBAL_RNG

"""Supertype of all symbolic operators. Subtype of `AbstractCliffordOperator`"""
abstract type AbstractSymbolicOperator <: AbstractCliffordOperator end
"""Supertype of all single-qubit symbolic operators."""
abstract type AbstractSingleQubitOperator <: AbstractSymbolicOperator end
"""Supertype of all two-qubit symbolic operators."""
abstract type AbstractTwoQubitOperator <: AbstractSymbolicOperator end
"""Supertype of all symbolic single-qubit measurements."""
abstract type AbstractMeasurement <: AbstractOperation end

const MINBATCH1Q = 100
const MINBATCH2Q = 100

# Stim has a good list of specialized single and two qubit operations at https://github.com/quantumlib/Stim/blob/e51ea66d213b25920e72c08e53266ec56fd14db4/src/stim/stabilizers/tableau_specialized_prepend.cc
# Note that their specialized operations are for prepends (right multiplications), while we implement append (left multiplication) operations.

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

##############################
# Single-qubit gates
##############################

function _apply!(stab::AbstractStabilizer, gate::G; phases::Val{B}=Val(true)) where {B, G<:AbstractSingleQubitOperator}
    s = tab(stab)
    c = gate.q
    @batch per=core minbatch=MINBATCH1Q for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        x,z,phase = qubit_kernel(gate,x,z)
        setxbit(s, r, c, x)
        setzbit(s, r, c, z)
        B && phase && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

"""Macro used to define single qubit symbolic gates and their `qubit_kernel` methods."""
macro qubitop1(name, kernel)
    prefixname = Symbol(:s,name)
    docstring = "A \"symbolic\" single-qubit $name. See also: [`SingleQubitOperator`](@ref), [`AbstractSymbolicOperator`](@ref)"
    quote
        struct $(esc(prefixname)) <: AbstractSingleQubitOperator
            q::Int
        end
        @doc $docstring $prefixname
        @inline $(esc(:qubit_kernel))(::$prefixname, x, z) = $kernel
    end
end

@qubitop1 Hadamard (z , x   , x!=0 && z!=0)
@qubitop1 Phase    (x , x⊻z , x!=0 && z!=0)
@qubitop1 InvPhase (x , x⊻z , x!=0 && z==0)
@qubitop1 X        (x , z   , z!=0)
@qubitop1 Y        (x , z   , (x⊻z)!=0)
@qubitop1 Z        (x , z   , x!=0)

"""A "symbolic" single-qubit Identity operation.

See also: [`SingleQubitOperator`](@ref)
"""
struct sId1 <: AbstractSingleQubitOperator
    q::Int
end
function _apply!(stab::AbstractStabilizer, ::sId1; phases::Val{B}=Val(true)) where B
    stab
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
CliffordOperator{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}}

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
function _apply!(stab::AbstractStabilizer, op::SingleQubitOperator; phases::Val{B}=Val(true)) where B # TODO Generated functions that simplify the whole `if phases` branch might be a good optimization, but a quick benchmakr comparing sHadamard to SingleQubitOperator(sHadamard) did not show a worthwhile difference.
    s = tab(stab)
    c = op.q
    Tme = eltype(s.xzs)
    sh = getshift(Tme, c)
    xx,zx,xz,zz = Tme.((op.xx,op.zx,op.xz,op.zz)) .<< sh
    anticom = ~iszero((~zz & xz & ~xx & zx) | ( zz & ~xz & xx & zx) | (zz &  xz & xx & ~zx))
    @batch per=core minbatch=MINBATCH1Q for r in eachindex(s)
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
        n==1 || throw(ArgumentError("Set `n=1` as a `SingleQubitOperator` being compacted (`compact=true`) has to result in a 1×1 `CliffordOperator`."))
        return CliffordOperator(Tableau([op.px ? 0x2 : 0x0, op.pz ? 0x2 : 0x0],[op.xx op.xz; op.zx op.zz]))
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

CliffordOperator(::Type{O}) where {O<:AbstractSingleQubitOperator} = CliffordOperator(apply!(one(Destabilizer,1),O(1)))

function Base.show(io::IO, op::AbstractSingleQubitOperator)
    print(io, "Symbolic single-qubit gate on qubit $(op.q)\n")
    show(io, CliffordOperator(op,1;compact=true))
end

"""Random symbolic single-qubit Clifford applied to qubit at index `qubit`.

See also: [`SingleQubitOperator`](@ref), [`random_clifford`](@ref)
"""
function random_clifford1(rng::AbstractRNG, qubit)
    return enumerate_single_qubit_gates(rand(rng,1:6),qubit=qubit,phases=(rand(rng,Bool),rand(rng,Bool)))
end
random_clifford1(qubit) = random_clifford1(GLOBAL_RNG, qubit)

##############################
# Two-qubit gates
##############################

function _apply!(stab::AbstractStabilizer, gate::G; phases::Val{B}=Val(true)) where {B, G<:AbstractTwoQubitOperator}
    s = tab(stab)
    q1 = gate.q1
    q2 = gate.q2
    Tme = eltype(s.xzs)
    shift = getshift(Tme, q1) - getshift(Tme, q2)
    @batch per=core minbatch=MINBATCH2Q for r in eachindex(s)
        x1 = getxbit(s, r, q1)
        z1 = getzbit(s, r, q1)
        x2 = getxbit(s, r, q2)<<shift
        z2 = getzbit(s, r, q2)<<shift
        x1,z1,x2,z2,phase = qubit_kernel(gate,x1,z1,x2,z2) # Most `qubit_kernel` functions are defined by a `qubitop2` macro
        setxbit(s, r, q1, x1, 0)
        setzbit(s, r, q1, z1, 0)
        setxbit(s, r, q2, x2, -shift)
        setzbit(s, r, q2, z2, -shift)
        if B && phase
            s.phases[r] += 0x2
            s.phases[r] &= 3
        end
    end
    stab
end

"""Macro used to define 2-qubit symbolic gates and their `qubit_kernel` methods."""
macro qubitop2(name, kernel)
    prefixname = Symbol(:s,name)
    docstring = "A \"symbolic\" $name. See also: [`AbstractSymbolicOperator`](@ref)"
    quote
        struct $(esc(prefixname)) <: AbstractTwoQubitOperator
            q1::Int
            q2::Int
        end
        @doc $docstring $prefixname
        @inline $(esc(:qubit_kernel))(::$prefixname, x1, z1, x2, z2) = $kernel
    end
end

@qubitop2 SWAP   (x2 , z2    , x1    , z1    , false)
@qubitop2 CNOT   (x1 , z1⊻z2 , x2⊻x1 , z2    , ~iszero( (x1 & z1 & x2 & z2)  | (x1 & z2 &~(z1|x2)) ))
@qubitop2 CPHASE (x1 , z1⊻x2 , x2    , z2⊻x1 , ~iszero( (x1 & z1 & x2 &~z2)  | (x1 &~z1 & x2 & z2) ))

function CliffordOperator(op::AbstractTwoQubitOperator, n; compact=false)
    if compact
        n==2 || throw(ArgumentError("Set `n=2` as a `TwoQubitOperator` being compacted (`compact=true`) has to result in a 2×2 `CliffordOperator`."))
        return CliffordOperator(typeof(op))
    else
        n >= max(op.q1,op.q2) || throw(DimensionMismatch("Set a larger `n`, otherwise the `TwoQubitOperator` can not fit in the allocated `CliffordOperator`. Use `compact=true` if you want to discard the target index."))
        c = one(CliffordOperator, n)
        _c = CliffordOperator(typeof(op))
        for (i,q) in ((1,op.q1),(2,op.q2))
            for (ii,qq) in ((1,op.q1),(2,op.q2))
                c[q,qq] = _c[i,ii] # TODO define an `embed` helper function
                c[n+q,qq] = _c[n+i,ii]
            end
            c.tab.phases[q] = _c.tab.phases[i] # TODO define a `phasesview` or `phases` helper function
            c.tab.phases[n+q] = _c.tab.phases[n+i]
        end
        return c
    end
end

CliffordOperator(::Type{O}) where {O<:AbstractTwoQubitOperator} = CliffordOperator(apply!(one(Destabilizer,2),O(1,2)))

function Base.show(io::IO, op::AbstractTwoQubitOperator)
    print(io, "Symbolic two-qubit gate on qubit $(op.q1) and $(op.q2)\n")
    show(io, CliffordOperator(typeof(op)))
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

##############################
# Measurements
##############################

"""Symbolic single qubit X measurement. See also [`Register`](@ref), [`projectXrand!`](@ref), [`sMY`](@ref), [`sMZ`](@ref)"""
struct sMX{T<:Union{Int,Nothing}} <: AbstractMeasurement
    qubit::Int
    bit::T
end

"""Symbolic single qubit Y measurement. See also [`Register`](@ref), [`projectYrand!`](@ref), [`sMX`](@ref), [`sMZ`](@ref)"""
struct sMY{T<:Union{Int,Nothing}} <: AbstractMeasurement
    qubit::Int
    bit::T
end

"""Symbolic single qubit Z measurement. See also [`Register`](@ref), [`projectZrand!`](@ref), [`sMX`](@ref), [`sMY`](@ref)"""
struct sMZ{T<:Union{Int,Nothing}} <: AbstractMeasurement
    qubit::Int
    bit::T
end

sMX(i) = sMX(i,nothing)
sMY(i) = sMY(i,nothing)
sMZ(i) = sMZ(i,nothing)

function apply!(state::AbstractStabilizer, m::sMX)
    projectXrand!(state,m.qubit)
    state
end
function apply!(state::AbstractStabilizer, m::sMY)
    projectYrand!(state,m.qubit)
    state
end
function apply!(state::AbstractStabilizer, m::sMZ)
    projectZrand!(state,m.qubit)
    state
end
project!(state::AbstractStabilizer, m::sMX) = projectX!(state, m.qubit)
project!(state::AbstractStabilizer, m::sMY) = projectY!(state, m.qubit)
project!(state::AbstractStabilizer, m::sMZ) = projectZ!(state, m.qubit)
projectrand!(state::AbstractStabilizer, m::sMX) = projectXrand!(state, m.qubit)
projectrand!(state::AbstractStabilizer, m::sMY) = projectYrand!(state, m.qubit)
projectrand!(state::AbstractStabilizer, m::sMZ) = projectZrand!(state, m.qubit)
