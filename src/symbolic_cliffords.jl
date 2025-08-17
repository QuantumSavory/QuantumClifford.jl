using Random: AbstractRNG, GLOBAL_RNG
include("throws.jl")

THROW_COMPACT_QUBIT_OPERATOR_SHAPE_ERROR(n, operator_with_size) =
"Set `n=$n` as a $operator_with_size being compacted (`compact=true`) has to result in a $n x $n `CliffordOperator`."

THROW_INVALID_QUBIT_OPERATOR_SIZE_ALLOCATION(operator_with_size) =
"Set a larger `n`, otherwise the `operator_with_size` can not fit in the allocated `CliffordOperator`.
 Use `compact=true` if you want to discard the target index."

"""Supertype of all symbolic operators. Subtype of `AbstractCliffordOperator`"""
abstract type AbstractSymbolicOperator <: AbstractCliffordOperator end
"""Supertype of all single-qubit symbolic operators."""
abstract type AbstractSingleQubitOperator <: AbstractSymbolicOperator end
"""Supertype of all two-qubit symbolic operators."""
abstract type AbstractTwoQubitOperator <: AbstractSymbolicOperator end
"""Supertype of all symbolic single-qubit measurements."""
abstract type AbstractMeasurement <: AbstractOperation end

# Stim has a good list of specialized single and two qubit operations at https://github.com/quantumlib/Stim/blob/e51ea66d213b25920e72c08e53266ec56fd14db4/src/stim/stabilizers/tableau_specialized_prepend.cc
# Note that their specialized operations are for prepends (right multiplications), while we implement append (left multiplication) operations.

@inline getshift(::Type{Tₘₑ},col::Int) where {Tₘₑ} = _mod(Tₘₑ,col-1)
@inline getmask(::Type{Tₘₑ},col::Int) where {Tₘₑ} = Tₘₑ(0x1)<<getshift(Tₘₑ,col)
@inline getbigindex(::Type{Tₘₑ},col::Int) where {Tₘₑ} = _div(Tₘₑ,col-1)+1

TableauType{Tₚᵥ, Tₘₑ} = Tableau{Tₚᵥ, Tₘ} where {Tₘ <: AbstractMatrix{Tₘₑ}}

Base.@propagate_inbounds function getxbit(s::TableauType{Tₚᵥ, Tₘₑ}, r::Int, c::Int) where {Tₚᵥ, Tₘₑ}
    getxbit(s.xzs, r, c)
end
Base.@propagate_inbounds function getzbit(s::TableauType{Tₚᵥ, Tₘₑ}, r::Int, c::Int) where {Tₚᵥ, Tₘₑ}
    getzbit(s.xzs, r, c)
end
Base.@propagate_inbounds function setxbit(s::TableauType{Tₚᵥ, Tₘₑ}, r::Int, c::Int, x::Tₘₑ) where {Tₚᵥ, Tₘₑ}
    setxbit(s.xzs, r, c, x)
end
Base.@propagate_inbounds function setzbit(s::TableauType{Tₚᵥ, Tₘₑ}, r::Int, c::Int, z::Tₘₑ) where {Tₚᵥ, Tₘₑ}
    setzbit(s.xzs, r, c, z)
end
Base.@propagate_inbounds setxbit(s::TableauType{Tₚᵥ, Tₘₑ}, r::Int, c::Int, x::Tₘₑ, shift::Int) where {Tₚᵥ, Tₘₑ} = setxbit(s, r, c, x<<shift)
Base.@propagate_inbounds setzbit(s::TableauType{Tₚᵥ, Tₘₑ}, r::Int, c::Int, z::Tₘₑ, shift::Int) where {Tₚᵥ, Tₘₑ} = setzbit(s, r, c, z<<shift)


Base.@propagate_inbounds function getxbit(xzs::AbstractMatrix{T}, r::Int, c::Int)::T where {T <: Unsigned}
    xzs[QuantumClifford.getbigindex(T, c),r] & QuantumClifford.getmask(T, c)
end
Base.@propagate_inbounds function getzbit(xzs::AbstractMatrix{T}, r::Int, c::Int)::T where {T <: Unsigned}
    xzs[end÷2+QuantumClifford.getbigindex(T, c),r]& QuantumClifford.getmask(T, c)
end
Base.@propagate_inbounds function setxbit(xzs::AbstractMatrix{T}, r::Int, c::Int, x::T) where {T <: Unsigned}
    cbig = QuantumClifford.getbigindex(T, c)
    xzs[cbig,r] &= ~QuantumClifford.getmask(T, c)
    xzs[cbig,r] |= x
end
Base.@propagate_inbounds function setzbit(xzs::AbstractMatrix{T}, r::Int, c::Int, z::T) where {T <: Unsigned}
    cbig = QuantumClifford.getbigindex(T, c)
    xzs[end÷2+cbig,r] &= ~QuantumClifford.getmask(T, c)
    xzs[end÷2+cbig,r] |= z
end
Base.@propagate_inbounds setxbit(xzs::AbstractMatrix{T}, r::Int, c::Int, x::T, shift::Int) where {T <: Unsigned} = setxbit(xzs, r, c, x<<shift)
Base.@propagate_inbounds setzbit(xzs::AbstractMatrix{T}, r::Int, c::Int, z::T, shift::Int) where {T <: Unsigned} = setzbit(xzs, r, c, z<<shift)

##############################
# Single-qubit gates
##############################

function _apply!(stab::AbstractStabilizer, gate::AbstractSingleQubitOperator; phases::Val{B}=Val(true)) where {B}
    s = tab(stab)
    c = gate.q
    @inbounds @simd for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        x,z,phase = qubit_kernel(gate,x,z)
        setxbit(s, r, c, x)
        setzbit(s, r, c, z)
        B && phase && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function _apply_inv!(stab::AbstractStabilizer, gate::AbstractSingleQubitOperator; phases::Val{B}=Val(true)) where {B} # code repetition with the corresponding `_apply`
    s = tab(stab)
    c = gate.q
    @inbounds @simd for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        x,z,phase = inv_qubit_kernel(gate,x,z)
        setxbit(s, r, c, x)
        setzbit(s, r, c, z)
        B && phase && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

"""Macro used to define single qubit symbolic gates and their `qubit_kernel` methods."""
macro qubitop1(name, kernel, inv_kernel)
    prefixname = Symbol(:s,name)
    docstring = "A \"symbolic\" single-qubit $name. See also: [`SingleQubitOperator`](@ref), [`AbstractSymbolicOperator`](@ref)"
    quote
        struct $(esc(prefixname)) <: AbstractSingleQubitOperator
            q::Int
            $(esc(prefixname))(q) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q) end
        end
        @doc $docstring $prefixname
        @inline $(esc(:qubit_kernel))(::$prefixname, x, z) = $kernel
        @inline $(esc(:inv_qubit_kernel))(::$prefixname, x, z) = $inv_kernel
    end
end

@qubitop1 Hadamard     (z   ,x   , x!=0 && z!=0) (z   ,x   , x!=0 && z!=0)
@qubitop1 HadamardXY   (x   ,x⊻z , x==0 && z!=0) (x   ,x⊻z , x==0 && z!=0)
@qubitop1 HadamardYZ   (x⊻z ,z   , x!=0 && z==0) (x⊻z ,z   , x!=0 && z==0)
@qubitop1 Phase        (x   ,x⊻z , x!=0 && z!=0) (x   ,x⊻z , x!=0 && z==0)
@qubitop1 InvPhase     (x   ,x⊻z , x!=0 && z==0) (x   ,x⊻z , x!=0 && z!=0)
@qubitop1 X            (x   ,z   , z!=0)         (x   ,z   , z!=0)
@qubitop1 Y            (x   ,z   , (x⊻z)!=0)     (x   ,z   , (x⊻z)!=0)
@qubitop1 Z            (x   ,z   , x!=0)         (x   ,z   , x!=0)
@qubitop1 SQRTX        (x⊻z ,z   , x==0 && z!=0) (x⊻z ,z   , x!=0 && z!=0)
@qubitop1 InvSQRTX     (x⊻z ,z   , x!=0 && z!=0) (x⊻z ,z   , x==0 && z!=0)
@qubitop1 SQRTY        (z   ,x   , x!=0 && z==0) (z   ,x   , z!=0 && x==0)
@qubitop1 InvSQRTY     (z   ,x   , z!=0 && x==0) (z   ,x   , x!=0 && z==0) 
@qubitop1 CXYZ         (x⊻z ,x   , false)        (z   ,x⊻z , false)
@qubitop1 CZYX         (z   ,x⊻z , false)        (x⊻z ,x   , false)

"""A "symbolic" single-qubit Identity operation.

See also: [`SingleQubitOperator`](@ref)
"""
struct sId1 <: AbstractSingleQubitOperator
    q::Int
    sId1(q) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q) end
end
function _apply!(stab::AbstractStabilizer, ::sId1; phases::Val{B}=Val(true)) where B
    stab
end
function _apply_inv!(stab::AbstractStabilizer, ::sId1; phases::Val{B}=Val(true)) where B
    stab
end

"""A "symbolic" general single-qubit operator which permits faster multiplication than an operator expressed as an explicit tableau.

```jldoctest
julia> op = SingleQubitOperator(2, true, true, true, false, true, true) # Tableau components and phases
SingleQubitOperator on qubit 2
X₁ ⟼ - Y
Z₁ ⟼ - X

julia> typeof(op)
SingleQubitOperator

julia> t_op = CliffordOperator(op, 3) # Transforming it back into an explicit tableau representation (specifying the size)
X₁ ⟼ + X__
X₂ ⟼ - _Y_
X₃ ⟼ + __X
Z₁ ⟼ + Z__
Z₂ ⟼ - _X_
Z₃ ⟼ + __Z

julia> typeof(t_op)
CliffordOperator{QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}, PauliOperator{Array{UInt8, 0}, Vector{UInt64}}}

julia> CliffordOperator(op, 1, compact=true) # You can also extract just the non-trivial part of the tableau
X₁ ⟼ - Y
Z₁ ⟼ - X
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
    SingleQubitOperator(q,args...) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q,args...) end
end
function _apply!(stab::AbstractStabilizer, op::SingleQubitOperator; phases::Val{B}=Val(true)) where B # TODO Generated functions that simplify the whole `if phases` branch might be a good optimization, but a quick benchmakr comparing sHadamard to SingleQubitOperator(sHadamard) did not show a worthwhile difference.
    s = tab(stab)
    c = op.q
    Tₘₑ = eltype(s.xzs)
    sh = getshift(Tₘₑ, c)
    xx,zx,xz,zz = Tₘₑ.((op.xx,op.zx,op.xz,op.zz)) .<< sh
    anticom = ~iszero((~zz & xz & ~xx & zx) | ( zz & ~xz & xx & zx) | (zz &  xz & xx & ~zx))
    @inbounds @simd for r in eachindex(s)
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

SingleQubitOperator(h::sHadamard)           = SingleQubitOperator(h.q, false, true , true , false, false, false)
SingleQubitOperator(p::sPhase)              = SingleQubitOperator(p.q, true , true , false, true , false, false)
SingleQubitOperator(p::sInvPhase)           = SingleQubitOperator(p.q, true , true , false, true , true , false)
SingleQubitOperator(p::sId1)                = SingleQubitOperator(p.q, true , false, false, true , false, false)
SingleQubitOperator(p::sX)                  = SingleQubitOperator(p.q, true , false, false, true , false, true)
SingleQubitOperator(p::sY)                  = SingleQubitOperator(p.q, true , false, false, true , true , true)
SingleQubitOperator(p::sZ)                  = SingleQubitOperator(p.q, true , false, false, true , true , false)
SingleQubitOperator(p::sCXYZ)               = SingleQubitOperator(p.q, true , true , true , false, false, false)
SingleQubitOperator(p::sCZYX)               = SingleQubitOperator(p.q, false, true , true , true , false, false)
SingleQubitOperator(p::sHadamardXY)         = SingleQubitOperator(p.q, true , true , false, true , false, true)
SingleQubitOperator(p::sHadamardYZ)         = SingleQubitOperator(p.q, true , false, true , true , true , false)
SingleQubitOperator(p::sSQRTX)              = SingleQubitOperator(p.q, true , false, true , true , false, true)
SingleQubitOperator(p::sInvSQRTX)           = SingleQubitOperator(p.q, true , false, true , true , false, false)
SingleQubitOperator(p::sSQRTY)              = SingleQubitOperator(p.q, false, true , true , false, true , false)
SingleQubitOperator(p::sInvSQRTY)           = SingleQubitOperator(p.q, false, true , true , false, false, true)
SingleQubitOperator(o::SingleQubitOperator) = o
function SingleQubitOperator(o::SingleQubitOperator, qubit::Int)
    return SingleQubitOperator(qubit, o.xx, o.xz, o.zx, o.zz, o.px, o.pz)
end
function SingleQubitOperator(op::CliffordOperator, qubit)
    nqubits(op)==1 || throw(DimensionMismatch("You are trying to convert a multiqubit `CliffordOperator` into a symbolic `SingleQubitOperator`."))
    SingleQubitOperator(qubit,tab(op)[1,1]...,tab(op)[2,1]...,(~).(iszero.(tab(op).phases))...)
end
SingleQubitOperator(op::CliffordOperator) = SingleQubitOperator(op, 1)

CliffordOperator(op::AbstractSingleQubitOperator, n; kw...) = CliffordOperator(SingleQubitOperator(op), n; kw...)
function CliffordOperator(op::SingleQubitOperator, n; compact=false)
    if compact
        n==1 || throw(ArgumentError(THROW_COMPACT_QUBIT_OPERATOR_SHAPE_ERROR(1, "SingleQubitOperator")))
        return CliffordOperator(Tableau([op.px ? 0x2 : 0x0, op.pz ? 0x2 : 0x0],[op.xx op.xz; op.zx op.zz]))
    else
        n >= op.q || throw(DimensionMismatch(THROW_INVALID_QUBIT_OPERATOR_SIZE_ALLOCATION("SingleQubitOperator")))
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
    if get(io, :compact, false) | haskey(io, :typeinfo)
        print(io, "$(string(typeof(op)))($(op.q))")
    else
        print(io, "$(string(typeof(op))) on qubit $(op.q)\n")
        show(io, CliffordOperator(op,1;compact=true))
    end
end

"""Random symbolic single-qubit Clifford applied to qubit at index `qubit`.

See also: [`SingleQubitOperator`](@ref), [`random_clifford`](@ref)
"""
function random_clifford1(rng::AbstractRNG, qubit)
    return enumerate_single_qubit_gates(rand(rng,1:6),qubit=qubit,phases=(rand(rng,Bool),rand(rng,Bool)))
end
random_clifford1(qubit) = random_clifford1(GLOBAL_RNG, qubit)

function LinearAlgebra.inv(op::SingleQubitOperator)
    c = LinearAlgebra.inv(CliffordOperator(SingleQubitOperator(op), 1, compact=true))
    return SingleQubitOperator(c, op.q)
end

LinearAlgebra.inv(h::sHadamard)   = sHadamard(h.q)
LinearAlgebra.inv(p::sPhase)      = sInvPhase(p.q)
LinearAlgebra.inv(p::sInvPhase)   = sPhase(p.q)
LinearAlgebra.inv(p::sId1)        = sId1(p.q)
LinearAlgebra.inv(p::sX)          = sX(p.q)
LinearAlgebra.inv(p::sY)          = sY(p.q)
LinearAlgebra.inv(p::sZ)          = sZ(p.q)
LinearAlgebra.inv(p::sHadamardXY) = sHadamardXY(p.q)
LinearAlgebra.inv(p::sHadamardYZ) = sHadamardYZ(p.q)
LinearAlgebra.inv(p::sSQRTX)      = sInvSQRTX(p.q)
LinearAlgebra.inv(p::sInvSQRTX)   = sSQRTX(p.q)
LinearAlgebra.inv(p::sSQRTY)      = sInvSQRTY(p.q)
LinearAlgebra.inv(p::sInvSQRTY)   = sSQRTY(p.q)
LinearAlgebra.inv(p::sCZYX)       = sCXYZ(p.q)
LinearAlgebra.inv(p::sCXYZ)       = sCZYX(p.q)
##############################
# Two-qubit gates
##############################

function _apply!(stab::AbstractStabilizer, gate::AbstractTwoQubitOperator; phases::Val{B}=Val(true)) where {B}
    s = tab(stab)
    q1 = gate.q1
    q2 = gate.q2
    Tₘₑ = eltype(s.xzs)
    shift = getshift(Tₘₑ, q1) - getshift(Tₘₑ, q2)
    @inbounds @simd for r in eachindex(s)
#    for r in eachindex(s)
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

function _apply_inv!(stab::AbstractStabilizer, gate::AbstractTwoQubitOperator; phases::Val{B}=Val(true)) where {B} # code repetition with the corresponding `_apply`
    s = tab(stab)
    q1 = gate.q1
    q2 = gate.q2
    Tₘₑ = eltype(s.xzs)
    shift = getshift(Tₘₑ, q1) - getshift(Tₘₑ, q2)
    @inbounds @simd for r in eachindex(s)
#    for r in eachindex(s)
        x1 = getxbit(s, r, q1)
        z1 = getzbit(s, r, q1)
        x2 = getxbit(s, r, q2)<<shift
        z2 = getzbit(s, r, q2)<<shift
        x1,z1,x2,z2,phase = inv_qubit_kernel(gate,x1,z1,x2,z2) # Most `inv_qubit_kernel` functions are defined by a `qubitop2` macro
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
macro qubitop2(name, kernel, inv_kernel)
    prefixname = Symbol(:s,name)
    docstring = "A \"symbolic\" $name. See also: [`AbstractSymbolicOperator`](@ref)"
    quote
        struct $(esc(prefixname)) <: AbstractTwoQubitOperator
            q1::Int
            q2::Int
            $(esc(prefixname))(q1,q2) = if q1<=0 || q2<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) elseif q1==q2 throw(ArgumentError("Failed to create a two qubit gate because the two qubits it acts on have the same index. The qubit indices have to be different.")) else new(q1,q2) end
        end
        @doc $docstring $prefixname
        @inline $(esc(:qubit_kernel))(::$prefixname, x1, z1, x2, z2) = $kernel
        @inline $(esc(:inv_qubit_kernel))(::$prefixname, x1, z1, x2, z2) = $inv_kernel
    end
end
#                 x1   z1      x2      z2
@qubitop2 SWAP   (x2 , z2    , x1    , z1    , false) (x2 , z2    , x1    , z1    , false)

@qubitop2 SWAPCX    (x2    , z2⊻z1 , x2⊻x1 , z1    , ~iszero((x1 & z1 & x2 & z2) | (~x1 & z1 & x2 & ~z2))) (x2⊻x1 , z2    , x1    , z2⊻z1 , ~iszero((x1 & z1 & x2 & z2) | ( x1 &~z1 &~x2 &  z2)))
@qubitop2 InvSWAPCX (x2⊻x1 , z2    , x1    , z2⊻z1 , ~iszero((x1 & z1 & x2 & z2) | ( x1 &~z1 &~x2 &  z2))) (x2    , z2⊻z1 , x2⊻x1 , z1    , ~iszero((x1 & z1 & x2 & z2) | (~x1 & z1 & x2 & ~z2)))

@qubitop2 ISWAP    (x2       , x1⊻z2⊻x2 , x1       , x1⊻x2⊻z1 , ~iszero((x1 & z1 & ~x2) | (~x1 & x2 & z2))) (x2       , x1⊻z2⊻x2 , x1       , x1⊻x2⊻z1 , ~iszero((x1 &~z1 & ~x2) | (~x1 & x2 &~z2)))
@qubitop2 InvISWAP (x2       , x1⊻z2⊻x2 , x1       , x1⊻x2⊻z1 , ~iszero((x1 &~z1 & ~x2) | (~x1 & x2 &~z2))) (x2       , x1⊻z2⊻x2 , x1       , x1⊻x2⊻z1 , ~iszero((x1 & z1 & ~x2) | (~x1 & x2 & z2)))

@qubitop2 CZSWAP (x2    , z2⊻x1 , x1    , x2⊻z1 , ~iszero((x1 & ~z1 & x2 & z2) | (x1 & z1 & x2 & ~z2))) (x2    , z2⊻x1 , x1    , x2⊻z1 , ~iszero((x1 & ~z1 & x2 & z2) | (x1 & z1 & x2 & ~z2)))
@qubitop2 CXSWAP (x2⊻x1 , z2    , x1    , z2⊻z1 , ~iszero((x1 & ~z1 &~x2 & z2) | (x1 & z1 & x2 &  z2))) (x2 , z2⊻z1    , x2⊻x1    , z1 , ~iszero((x1 &  z1 & x2 & z2) | (~x1 & z1 & x2 & ~z2)))

@qubitop2 CNOT   (x1 , z1⊻z2 , x2⊻x1 , z2    , ~iszero( (x1 & z1 & x2 & z2)  | (x1 & z2 &~(z1|x2)) )) (x1 , z1⊻z2 , x2⊻x1 , z2    , ~iszero( (x1 & z1 & x2 & z2)  | (x1 & z2 &~(z1|x2)) ))
@qubitop2 CPHASE (x1 , z1⊻x2 , x2    , z2⊻x1 , ~iszero( (x1 & z1 & x2 &~z2)  | (x1 &~z1 & x2 & z2) )) (x1 , z1⊻x2 , x2    , z2⊻x1 , ~iszero( (x1 & z1 & x2 &~z2)  | (x1 &~z1 & x2 & z2) ))

@qubitop2 ZCX    (x1      , z1⊻z2    , x2⊻x1 , z2      , ~iszero( ((x1 & z2) &~(z1 ⊻ x2)) )) (x1      , z1⊻z2    , x2⊻x1 , z2      , ~iszero( ((x1 & z2) &~(z1 ⊻ x2)) )) # equiv of CNOT[1, 2]
@qubitop2 ZCY    (x1      , x2⊻z1⊻z2 , x2⊻x1 , z2⊻x1   , ~iszero( (x1 & (x2 ⊻ z1) & (x2 ⊻ z2)) )) (x1      , x2⊻z1⊻z2 , x2⊻x1 , z2⊻x1   , ~iszero( (x1 & (x2 ⊻ z1) & (x2 ⊻ z2)) ))
@qubitop2 ZCZ    (x1      , z1⊻x2    , x2    , z2⊻x1   , ~iszero( ((z1 ⊻ z2) & (x1 & x2)) )) (x1      , z1⊻x2    , x2    , z2⊻x1   , ~iszero( ((z1 ⊻ z2) & (x1 & x2)) ))

@qubitop2 XCX    (z2⊻x1   , z1       , x2⊻z1 , z2      , ~iszero( (z1 & z2) & (x1 ⊻ x2) )) (z2⊻x1   , z1       , x2⊻z1 , z2      , ~iszero( (z1 & z2) & (x1 ⊻ x2) ))
@qubitop2 XCY    (x1⊻x2⊻z2, z1       , z1⊻x2 , z2⊻z1   , ~iszero( (z1 & (x2 ⊻ z2) & ~(x2 ⊻ x1)) )) (x1⊻x2⊻z2, z1       , z1⊻x2 , z2⊻z1   , ~iszero( (z1 & (x2 ⊻ z2) & ~(x2 ⊻ x1)) ))
@qubitop2 XCZ    (x1⊻x2   , z1       , x2    , z2⊻z1   , ~iszero( (x2 & z1) & ~(x1 ⊻ z2) )) (x1⊻x2   , z1       , x2    , z2⊻z1   , ~iszero( (x2 & z1) & ~(x1 ⊻ z2) )) # equiv to CNOT[2, 1]

@qubitop2 YCX    (x1⊻z2   , z2⊻z1    , x1⊻z1⊻x2 , z2      , ~iszero( (z2 & (x1 ⊻ z1) & ~(x2 ⊻ x1)) )) (x1⊻z2   , z2⊻z1    , x1⊻z1⊻x2 , z2      , ~iszero( (z2 & (x1 ⊻ z1) & ~(x2 ⊻ x1)) ))
@qubitop2 YCY    (x1⊻z2⊻x2, z1⊻x2⊻z2 , x1⊻x2⊻z1 , x1⊻z1⊻z2, ~iszero( (x1 & ~z1 & ~x2 & z2) | (~x1 & z1 & x2 & ~z2))) (x1⊻z2⊻x2, z1⊻x2⊻z2 , x1⊻x2⊻z1 , x1⊻z1⊻z2, ~iszero( (x1 & ~z1 & ~x2 & z2) | (~x1 & z1 & x2 & ~z2)))
@qubitop2 YCZ    (x1⊻x2   , x2⊻z1    , x2       , z2⊻x1⊻z1, ~iszero( (x2 & (x1 ⊻ z1) & (z2 ⊻ x1)) )) (x1⊻x2   , x2⊻z1    , x2       , z2⊻x1⊻z1, ~iszero( (x2 & (x1 ⊻ z1) & (z2 ⊻ x1)) ))

@qubitop2 ZCrY    (x1, x1⊻z1⊻x2⊻z2, x1⊻x2, x1⊻z2, ~iszero((x1 &~z1 & x2) | (x1 & ~z1 & ~z2) | (x1 & x2 & ~z2))) (x1, x1⊻z1⊻x2⊻z2, x1⊻x2, x1⊻z2, ~iszero((x1 & z1 &~x2) | (x1 &  z1 &  z2) | (x1 &~x2 &  z2)))
@qubitop2 InvZCrY (x1, x1⊻z1⊻x2⊻z2, x1⊻x2, x1⊻z2, ~iszero((x1 & z1 &~x2) | (x1 &  z1 &  z2) | (x1 &~x2 &  z2))) (x1, x1⊻z1⊻x2⊻z2, x1⊻x2, x1⊻z2, ~iszero((x1 &~z1 & x2) | (x1 & ~z1 & ~z2) | (x1 & x2 & ~z2)))

@qubitop2 SQRTZZ    (x1       , x1⊻x2⊻z1 , x2       , x1⊻z2⊻x2 , ~iszero((x1 & z1 & ~x2) | (~x1 & x2 & z2))) (x1       , x1⊻x2⊻z1 , x2       , x1⊻z2⊻x2 , ~iszero((x1 &~z1 & ~x2) | (~x1 & x2 &~z2)))
@qubitop2 InvSQRTZZ (x1       , x1⊻x2⊻z1 , x2       , x1⊻z2⊻x2 , ~iszero((x1 &~z1 & ~x2) | (~x1 & x2 &~z2))) (x1       , x1⊻x2⊻z1 , x2       , x1⊻z2⊻x2 , ~iszero((x1 & z1 & ~x2) | (~x1 & x2 & z2)))

@qubitop2 SQRTXX    (z1⊻z2⊻x1, z1      , z1⊻x2⊻z2, z2      , ~iszero((~x1 & z1 &~z2) | (~z1 &~x2 & z2))) (z1⊻z2⊻x1, z1      , z1⊻x2⊻z2, z2      , ~iszero(( x1 & z1 &~z2) | (~z1 & x2 & z2)))
@qubitop2 InvSQRTXX (z1⊻z2⊻x1, z1      , z1⊻x2⊻z2, z2      , ~iszero(( x1 & z1 &~z2) | (~z1 & x2 & z2))) (z1⊻z2⊻x1, z1      , z1⊻x2⊻z2, z2      , ~iszero((~x1 & z1 &~z2) | (~z1 &~x2 & z2)))

@qubitop2 SQRTYY    (z1⊻x2⊻z2, x1⊻z2⊻x2, x1⊻z1⊻z2, x1⊻x2⊻z1, ~iszero((~x1 &~z1 & x2 &~z2) | ( x1 &~z1 &~x2 &~z2) | ( x1 &~z1 & x2 & z2) | ( x1 & z1 & x2 &~z2))) (z1⊻x2⊻z2, x1⊻z2⊻x2, x1⊻z1⊻z2, x1⊻x2⊻z1, ~iszero(( x1 & z1 &~x2 & z2) | (~x1 & z1 & x2 & z2) | (~x1 & z1 &~x2 &~z2) | (~x1 &~z1 &~x2 & z2)))
@qubitop2 InvSQRTYY (z1⊻x2⊻z2, x1⊻z2⊻x2, x1⊻z1⊻z2, x1⊻x2⊻z1, ~iszero(( x1 & z1 &~x2 & z2) | (~x1 & z1 & x2 & z2) | (~x1 & z1 &~x2 &~z2) | (~x1 &~z1 &~x2 & z2))) (z1⊻x2⊻z2, x1⊻z2⊻x2, x1⊻z1⊻z2, x1⊻x2⊻z1, ~iszero((~x1 &~z1 & x2 &~z2) | ( x1 &~z1 &~x2 &~z2) | ( x1 &~z1 & x2 & z2) | ( x1 & z1 & x2 &~z2)))

#=
To get the boolean formulas for the phase, it is easiest to first write down the truth table for the phase:
for i in 0:15
    s=S"__"
    d=tuple(Bool.(digits(i,base=2,pad=4))...)
    s[1,1]=d[1:2]
    s[1,2]=d[3:4]
    println("$(Int.(d)) $(phases(explicit_tableaux*s)[1]÷2)")
end
Then we can use a truth-table to boolean formula converter, e.g. https://www.dcode.fr/boolean-expressions-calculator
(by just writing the initial unsimplified formula as the OR of all true rows of the table)
=#

function CliffordOperator(op::AbstractTwoQubitOperator, n; compact=false)
    if compact
        n==2 || throw(ArgumentError(THROW_COMPACT_QUBIT_OPERATOR_SHAPE_ERROR(2, "TwoQubitOperator")))
        return CliffordOperator(typeof(op))
    else
        n >= max(op.q1,op.q2) || throw(DimensionMismatch(THROW_INVALID_QUBIT_OPERATOR_SIZE_ALLOCATION("TwoQubitOperator")))
        c = one(CliffordOperator, n)
        _c = CliffordOperator(typeof(op))
        for (i,q) in ((1,op.q1),(2,op.q2))
            for (ii,qq) in ((1,op.q1),(2,op.q2))
                c[q,qq] = _c[i,ii] # TODO define an `embed` helper function
                c[n+q,qq] = _c[2+i,ii]
            end
            c.tab.phases[q] = _c.tab.phases[i] # TODO define a `phasesview` or `phases` helper function
            c.tab.phases[n+q] = _c.tab.phases[2+i]
        end
        return c
    end
end

CliffordOperator(::Type{O}) where {O<:AbstractTwoQubitOperator} = CliffordOperator(apply!(one(Destabilizer,2),O(1,2)))

function Base.show(io::IO, op::AbstractTwoQubitOperator)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        print(io, "$(string(typeof(op)))($(op.q1),$(op.q2))")
    else
        print(io, "$(string(typeof(op))) on qubit ($(op.q1),$(op.q2))\n")
        show(io, CliffordOperator(op,2;compact=true))
    end
end

LinearAlgebra.inv(op::sSWAP)      = sSWAP(op.q1, op.q2)
LinearAlgebra.inv(op::sCNOT)      = sCNOT(op.q1, op.q2)
LinearAlgebra.inv(op::sCPHASE)    = sCPHASE(op.q1, op.q2)
LinearAlgebra.inv(op::sZCX)       = sZCX(op.q1, op.q2)
LinearAlgebra.inv(op::sZCY)       = sZCY(op.q1, op.q2)
LinearAlgebra.inv(op::sZCZ)       = sZCZ(op.q1, op.q2)
LinearAlgebra.inv(op::sXCX)       = sXCX(op.q1, op.q2)
LinearAlgebra.inv(op::sXCY)       = sXCY(op.q1, op.q2)
LinearAlgebra.inv(op::sXCZ)       = sXCZ(op.q1, op.q2)
LinearAlgebra.inv(op::sYCX)       = sYCX(op.q1, op.q2)
LinearAlgebra.inv(op::sYCY)       = sYCY(op.q1, op.q2)
LinearAlgebra.inv(op::sYCZ)       = sYCZ(op.q1, op.q2)
LinearAlgebra.inv(op::sZCrY)      = sInvZCrY(op.q1, op.q2)
LinearAlgebra.inv(op::sInvZCrY)   = sZCrY(op.q1, op.q2)
LinearAlgebra.inv(op::sSWAPCX)    = sInvSWAPCX(op.q1, op.q2)
LinearAlgebra.inv(op::sInvSWAPCX) = sSWAPCX(op.q1, op.q2)
LinearAlgebra.inv(op::sCZSWAP)    = sCZSWAP(op.q1, op.q2)
LinearAlgebra.inv(op::sCXSWAP)    = sSWAPCX(op.q1, op.q2)
LinearAlgebra.inv(op::sISWAP)     = sInvISWAP(op.q1, op.q2)
LinearAlgebra.inv(op::sInvISWAP)  = sISWAP(op.q1, op.q2)
LinearAlgebra.inv(op::sSQRTZZ)    = sInvSQRTZZ(op.q1, op.q2)
LinearAlgebra.inv(op::sInvSQRTZZ) = sSQRTZZ(op.q1, op.q2)
LinearAlgebra.inv(op::sSQRTXX)    = sInvSQRTXX(op.q1, op.q2)
LinearAlgebra.inv(op::sInvSQRTXX) = sSQRTXX(op.q1, op.q2)
LinearAlgebra.inv(op::sSQRTYY)    = sInvSQRTYY(op.q1, op.q2)
LinearAlgebra.inv(op::sInvSQRTYY) = sSQRTYY(op.q1, op.q2)

##############################
# Functions that perform direct application of common operators without needing an operator instance
##############################
# TODO is there a reason to keep these given that sZ/sX/sY exist?
# TODO currently these are faster than sZ/sX/sY

"""Apply a Pauli Z to the `i`-th qubit of state `s`. You should use `apply!(stab,sZ(i))` instead of this."""
function apply_single_z!(stab::AbstractStabilizer, i)
    s = tab(stab)
    _, ibig, _, ismallm = get_bitmask_idxs(s.xzs,i)
    @inbounds @simd for row in 1:size(s.xzs,2)
        if !iszero(s.xzs[ibig,row] & ismallm)
            s.phases[row] = (s.phases[row]+0x2)&0x3
        end
    end
    stab
end

"""Apply a Pauli X to the `i`-th qubit of state `s`. You should use `apply!(stab,sX(i))` instead of this."""
function apply_single_x!(stab::AbstractStabilizer, i)
    s = tab(stab)
    _, ibig, _, ismallm = get_bitmask_idxs(s.xzs,i)
    @inbounds @simd for row in 1:size(s.xzs,2)
        if !iszero(s.xzs[end÷2+ibig,row] & ismallm)
            s.phases[row] = (s.phases[row]+0x2)&0x3
        end
    end
    stab
end

"""Apply a Pauli Y to the `i`-th qubit of state `s`. You should use `apply!(stab,sY(i))` instead of this."""
function apply_single_y!(stab::AbstractStabilizer, i)
    s = tab(stab)
    _, ibig, _, ismallm = get_bitmask_idxs(s.xzs,i)
    @inbounds @simd for row in 1:size(s.xzs,2)
        if !iszero((s.xzs[ibig,row] & ismallm) ⊻ (s.xzs[end÷2+ibig,row] & ismallm))
            s.phases[row] = (s.phases[row]+0x2)&0x3
        end
    end
    stab
end

##############################
# Measurements
##############################

"""Symbolic single qubit X measurement. See also [`Register`](@ref), [`projectXrand!`](@ref), [`sMY`](@ref), [`sMZ`](@ref)"""
struct sMX <: AbstractMeasurement
    qubit::Int
    bit::Int
    sMX(q, args...) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q,args...) end
end

"""Symbolic single qubit Y measurement. See also [`Register`](@ref), [`projectYrand!`](@ref), [`sMX`](@ref), [`sMZ`](@ref)"""
struct sMY <: AbstractMeasurement
    qubit::Int
    bit::Int
    sMY(q, args...) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q,args...) end
end

"""Symbolic single qubit Z measurement. See also [`Register`](@ref), [`projectZrand!`](@ref), [`sMX`](@ref), [`sMY`](@ref)"""
struct sMZ <: AbstractMeasurement
    qubit::Int
    bit::Int
    sMZ(q, args...) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q,args...) end
end

sMX(i) = sMX(i,0)
sMY(i) = sMY(i,0)
sMZ(i) = sMZ(i,0)
sMX(i,::Nothing) = sMX(i,0)
sMY(i,::Nothing) = sMY(i,0)
sMZ(i,::Nothing) = sMZ(i,0)

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

"""Measure a qubit in the Z basis and reset to the |0⟩ state.

!!! warning "It does not trace out the qubit!"
    As described below there is a difference between measuring the qubit (followed by setting it to a given known state)
    and "tracing out" the qubit. By reset here we mean "measuring and setting to a known state", not "tracing out".

```jldoctest
julia> s = MixedDestabilizer(S"XXX ZZI IZZ") # |000⟩+|111⟩
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z__
+ _X_
+ __X
𝒮𝓉𝒶𝒷━
+ XXX
+ ZZ_
+ Z_Z

julia> traceout!(copy(s), 1) # = I⊗(|00⟩⟨00| + |11⟩⟨11|)
𝒟ℯ𝓈𝓉𝒶𝒷
+ _X_
𝒳ₗ━━━
+ _XX
+ Z__
𝒮𝓉𝒶𝒷━
+ _ZZ
𝒵ₗ━━━
+ Z_Z
+ XXX

julia> projectZ!(traceout!(copy(s), 1), 1)[1] # = |000⟩⟨000|+|011⟩⟨011| or |100⟩⟨100|+|111⟩⟨111| (use projectZrand! to actually get a random result)
𝒟ℯ𝓈𝓉𝒶𝒷
+ _X_
+ XXX
𝒳ₗ━━━
+ _XX
𝒮𝓉𝒶𝒷━
+ _ZZ
+ Z__
𝒵ₗ━━━
+ Z_Z

julia> projectZ!(copy(s), 1)[1] # = |000⟩ or |111⟩ (use projectZrand! to actually get a random result)
𝒟ℯ𝓈𝓉𝒶𝒷
+ XXX
+ _X_
+ __X
𝒮𝓉𝒶𝒷━
+ Z__
+ ZZ_
+ Z_Z
```

```julia-repl
julia> apply!(Register(copy(s)), sMRZ(1)) |> quantumstate # |000⟩ or |011⟩, depending on randomization
𝒟ℯ𝓈𝓉𝒶𝒷
+ XXX
+ _X_
+ __X
𝒮𝓉𝒶𝒷━
+ Z__
- ZZ_
- Z_Z
```

See also: [`Reset`](@ref), [`sMZ`](@ref)"""
struct sMRZ <: AbstractOperation
    qubit::Int
    bit::Int
    sMRZ(q, args...) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q,args...) end
end

"""Measure a qubit in the X basis and reset to the |+⟩ state.

See also: [`sMRZ`](@ref), [`Reset`](@ref), [`sMZ`](@ref)"""
struct sMRX <: AbstractOperation
    qubit::Int
    bit::Int
    sMRX(q, args...) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q,args...) end
end

"""Measure a qubit in the Y basis and reset to the |i₊⟩ state.

See also: [`sMRZ`](@ref), [`Reset`](@ref), [`sMZ`](@ref)"""
struct sMRY <: AbstractOperation
    qubit::Int
    bit::Int
    sMRY(q, args...) = if q<=0 throw(ArgumentError(THROW_NO_ZERO_QUBIT)) else new(q,args...) end
end

sMRX(i) = sMRX(i,0)
sMRY(i) = sMRY(i,0)
sMRZ(i) = sMRZ(i,0)
sMRX(i,::Nothing) = sMRX(i,0)
sMRY(i,::Nothing) = sMRY(i,0)
sMRZ(i,::Nothing) = sMRZ(i,0)
