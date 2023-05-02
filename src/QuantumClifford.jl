"""
A module for using the Stabilizer formalism and simulating Clifford circuits.
"""
module QuantumClifford

# TODO document phases=false

# TODO Significant performance improvements: many operations do not need phase=true if the Pauli operations commute

import LinearAlgebra
using LinearAlgebra: inv, mul!, rank
using DocStringExtensions
using Polyester
#using LoopVectorization
using HostCPUFeatures: pick_vector_width
import SIMD

import QuantumInterface: tensor, ⊗, tensor_pow, apply!, nqubits, expect, project!, reset_qubits!, traceout!, ptrace, apply!, projectX!, projectY!, projectZ!, entanglement_entropy

export
    @P_str, PauliOperator, ⊗, I, X, Y, Z,
    @T_str, xbit, zbit, xview, zview,
    @S_str, Stabilizer,
    Destabilizer, MixedStabilizer, MixedDestabilizer,
    prodphase, comm,
    nqubits,
    stabilizerview, destabilizerview, logicalxview, logicalzview, phases,
    bitview, quantumstate, tab,
    BadDataStructure,
    affectedqubits,
    # GF2
    stab_to_gf2, gf2_gausselim!, gf2_isinvertible, gf2_invert, gf2_H_to_G,
    # Canonicalization
    canonicalize!, canonicalize_rref!, canonicalize_gott!,
    # Linear Algebra
    tensor, tensor_pow,
    logdot, expect,
    apply!,
    # Low Level Function Interface
    generate!, project!, reset_qubits!, traceout!,
    projectX!, projectY!, projectZ!,
    projectrand!, projectXrand!, projectYrand!, projectZrand!,
    puttableau!,
    # Clifford Ops
    CliffordOperator, @C_str, permute,
    tCNOT, tCPHASE, tSWAP, tHadamard, tPhase, tId1,
    # Symbolic Clifford Ops
    AbstractSymbolicOperator, AbstractSingleQubitOperator, AbstractTwoQubitOperator,
    sHadamard, sPhase, sInvPhase, SingleQubitOperator, sId1, sX, sY, sZ,
    sCNOT, sCPHASE, sSWAP,
    sXCX, sXCY, sXCZ, sYCX, sYCY, sYCZ, sZCX, sZCY, sZCZ,
    # Misc Ops
    SparseGate,
    sMX, sMY, sMZ, PauliMeasurement, Reset,
    BellMeasurement,
    Register,
    # Enumeration and Randoms
    enumerate_single_qubit_gates, random_clifford1,
    enumerate_cliffords, symplecticGS, clifford_cardinality, enumerate_phases,
    random_invertible_gf2,
    random_pauli, random_stabilizer, random_destabilizer, random_clifford,
    # Noise
    applynoise!, UnbiasedUncorrelatedNoise, NoiseOp, NoiseOpAll, NoisyGate,
    # Useful States
    bell, ghz,
    single_z, single_x, single_y,
    # Graphs
    graphstate, graphstate!, graph_gatesequence, graph_gate,
    # Clipped Gauge
    canonicalize_clip!, bigram, entanglement_entropy,
    # mctrajectories
    CircuitStatus, continue_stat, true_success_stat, false_success_stat, failure_stat,
    mctrajectory!, mctrajectories, applywstatus!


const BIG_INT_MINUS_ONE = Ref{BigInt}()
const BIG_INT_TWO = Ref{BigInt}()
const BIG_INT_FOUR = Ref{BigInt}()

function __init__()
    BIG_INT_MINUS_ONE[] = BigInt(-1)
    BIG_INT_TWO[] = BigInt(2)
    BIG_INT_FOUR[] = BigInt(4)
end

# Predefined constants representing the permitted phases encoded
# in the low bits of UInt8.
const _p  = 0x00
const _pi = 0x01
const _m  = 0x02
const _mi = 0x03

const phasedict = Dict(""=>_p,"+"=>_p,"i"=>_pi,"+i"=>_pi,"-"=>_m,"-i"=>_mi)
const toletter = Dict((false,false)=>"_",(true,false)=>"X",(false,true)=>"Z",(true,true)=>"Y")

include("macrotools.jl")

##############################
# Pauli Operators
##############################

abstract type AbstractOperation end
abstract type AbstractCliffordOperator <: AbstractOperation end

"""
A multi-qubit Pauli operator (``±\\{1,i\\}\\{I,Z,X,Y\\}^{\\otimes n}``).

A Pauli can be constructed with the `P` custom string macro or by building
up one through products and tensor products of smaller operators.

```jldoctest
julia> pauli3 = P"-iXYZ"
-iXYZ

julia> pauli4 = 1im * pauli3 ⊗ X
+ XYZX

julia> Z*X
+iY
```

We use a typical F(2,2) encoding internally. The X and Z bits are stored
in a single concatenated padded array of UInt chunks of a bit array.

```jldoctest
julia> p = P"-IZXY";


julia> p.xz
2-element Vector{UInt64}:
 0x000000000000000c
 0x000000000000000a
```

You can access the X and Z bits through getters and setters or through the
`xview`, `zview`, `xbit`, and `zbit` functions.

```jldoctest
julia> p = P"XYZ"; p[1]
(true, false)

julia> p[1] = (true, true); p
+ YYZ
```
"""
struct PauliOperator{Tz<:AbstractArray{UInt8,0}, Tv<:AbstractVector{<:Unsigned}} <: AbstractCliffordOperator
    phase::Tz
    nqubits::Int
    xz::Tv
end

PauliOperator(phase::UInt8, nqubits::Int, xz::Tv) where Tv<:AbstractVector{<:Unsigned} = PauliOperator(fill(UInt8(phase),()), nqubits, xz)
PauliOperator(phase::UInt8, x::AbstractVector{Bool}, z::AbstractVector{Bool}) = PauliOperator(fill(UInt8(phase),()), length(x), vcat(reinterpret(UInt,BitVector(x).chunks),reinterpret(UInt,BitVector(z).chunks)))
PauliOperator(x::AbstractVector{Bool}, z::AbstractVector{Bool}) = PauliOperator(0x0, x, z)
PauliOperator(xz::AbstractVector{Bool}) = PauliOperator(0x0, (@view xz[1:end÷2]), (@view xz[end÷2+1:end]))

"""Get a view of the X part of the `UInt` array of packed qubits of a given Pauli operator."""
function xview(p::PauliOperator)
    @view p.xz[1:end÷2]
end
"""Get a view of the Y part of the `UInt` array of packed qubits of a given Pauli operator."""
function zview(p::PauliOperator)
    @view p.xz[end÷2+1:end]
end
"""Extract as a new bit array the X part of the `UInt` array of packed qubits of a given Pauli operator."""
function xbit(p::PauliOperator)
    one = eltype(p.xz)(1)
    size = sizeof(eltype(p.xz))*8
    [(word>>s)&one==one for word in xview(p) for s in 0:size-1][begin:p.nqubits]
end
"""Extract as a new bit array the Z part of the `UInt` array of packed qubits of a given Pauli operator."""
function zbit(p::PauliOperator)
    one = eltype(p.xz)(1)
    size = sizeof(eltype(p.xz))*8
    [(word>>s)&one==one for word in zview(p) for s in 0:size-1][begin:p.nqubits]
end

function _P_str(a)
    letters = filter(x->occursin(x,"_IZXY"),a)
    phase = phasedict[strip(filter(x->!occursin(x,"_IZXY"),a))]
    PauliOperator(phase, [l=='X'||l=='Y' for l in letters], [l=='Z'||l=='Y' for l in letters])
end

macro P_str(a)
    quote _P_str($a) end
end

Base.getindex(p::PauliOperator{Tz,Tv}, i::Int) where {Tz, Tve<:Unsigned, Tv<:AbstractVector{Tve}} = (p.xz[_div(Tve, i-1)+1] & Tve(0x1)<<_mod(Tve,i-1))!=0x0, (p.xz[end÷2+_div(Tve,i-1)+1] & Tve(0x1)<<_mod(Tve,i-1))!=0x0
Base.getindex(p::PauliOperator{Tz,Tv}, r) where {Tz, Tve<:Unsigned, Tv<:AbstractVector{Tve}} = PauliOperator(p.phase[], xbit(p)[r], zbit(p)[r])

function Base.setindex!(p::PauliOperator{Tz,Tv}, (x,z)::Tuple{Bool,Bool}, i) where {Tz, Tve, Tv<:AbstractVector{Tve}}
    if x
        p.xz[_div(Tve,i-1)+1] |= Tve(0x1)<<_mod(Tve,i-1)
    else
        p.xz[_div(Tve,i-1)+1] &= ~(Tve(0x1)<<_mod(Tve,i-1))
    end
    if z
        p.xz[end÷2+_div(Tve,i-1)+1] |= Tve(0x1)<<_mod(Tve,i-1)
    else
        p.xz[end÷2+_div(Tve,i-1)+1] &= ~(Tve(0x1)<<_mod(Tve,i-1))
    end
    p
end

Base.firstindex(p::PauliOperator) = 1

Base.lastindex(p::PauliOperator) = p.nqubits

Base.eachindex(p::PauliOperator) = 1:p.nqubits

Base.size(pauli::PauliOperator) = (pauli.nqubits,)

Base.length(pauli::PauliOperator) = pauli.nqubits

nqubits(pauli::PauliOperator) = pauli.nqubits

xz2str(x,z) = join(toletter[e] for e in zip(x,z))

function xz2str_limited(x,z, limit=50)
    tupl = collect(zip(x,z))
    n = length(tupl)
    if (limit >= n || limit == -1)
        return xz2str(x, z)
    end
    padding = limit÷2
    return join(toletter[tupl[i]] for i in 1:padding) * "…" * join(toletter[tupl[i]] for i in (n-padding):n)
end

_show(io::IO, p::PauliOperator, limit=50) = print(io, ["+ ","+i","- ","-i"][p.phase[]+1]*xz2str_limited(xbit(p),zbit(p), limit))

function Base.show(io::IO, p::PauliOperator)
    if get(io, :compact, false)
        _show(io, p, 10)
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, p, sz[2]-7)
    else
        _show(io, p, -1)
    end
end

Base.:(==)(l::PauliOperator, r::PauliOperator) = r.phase==l.phase && r.nqubits==l.nqubits && r.xz==l.xz

Base.hash(p::PauliOperator, h::UInt) = hash(p.phase,hash(p.nqubits,hash(p.xz, h)))

Base.copy(p::PauliOperator) = PauliOperator(copy(p.phase),p.nqubits,copy(p.xz))

_nchunks(i::Int,T::Type{<:Unsigned}) = 2*( (i-1) ÷ (8*sizeof(T)) + 1 )
Base.zero(::Type{PauliOperator{Tz, Tv}}, q) where {Tz,T<:Unsigned,Tv<:AbstractVector{T}} = PauliOperator(zeros(UInt8), q, zeros(T, _nchunks(q,T)))
Base.zero(::Type{PauliOperator}, q) = zero(PauliOperator{Array{UInt8, 0}, Vector{UInt}}, q)
Base.zero(p::P) where {P<:PauliOperator} = zero(P, nqubits(p))

"""Zero-out the phases and single-qubit operators in a [`PauliOperator`](@ref)"""
@inline function zero!(p::PauliOperator{Tz,Tv}) where {Tz, Tve<:Unsigned, Tv<:AbstractVector{Tve}}
    fill!(p.xz, zero(Tve))
    p.phase[] = 0x0
    p
end

##############################
# Generic Tableaux
##############################

struct Tableau{Tzv<:AbstractVector{UInt8}, Tm<:AbstractMatrix{<:Unsigned}}
    phases::Tzv
    nqubits::Int
    xzs::Tm
end

function Tableau(paulis::AbstractVector{PauliOperator{Tz,Tv}}) where {Tz<:AbstractArray{UInt8,0},Tve<:Unsigned,Tv<:AbstractVector{Tve}}
    r = length(paulis)
    n = nqubits(paulis[1])
    tab = zero(Tableau{Vector{UInt8},Matrix{Tve}},r,n)
    for i in eachindex(paulis)
        tab[i] = paulis[i]
    end
    tab
end

Tableau(phases::AbstractVector{UInt8}, xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Tableau(
    phases, size(xs,2),
    vcat(hcat((BitArray(xs[i,:]).chunks for i in 1:size(xs,1))...)::Matrix{UInt},
         hcat((BitArray(zs[i,:]).chunks for i in 1:size(zs,1))...)::Matrix{UInt}) # type assertions to help Julia infer types
)

Tableau(phases::AbstractVector{UInt8}, xzs::AbstractMatrix{Bool}) = Tableau(phases, xzs[:,1:end÷2], xzs[:,end÷2+1:end])

Tableau(xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Tableau(zeros(UInt8, size(xs,1)), xs, zs)

Tableau(xzs::AbstractMatrix{Bool}) = Tableau(zeros(UInt8, size(xzs,1)), xzs[:,1:end÷2], xzs[:,end÷2+1:end])

Tableau(t::Tableau) = t

function _T_str(a) # TODO this can be optimized by not creating intermediary PauliOperator objects
    paulis = [_P_str(strip(s.match)) for s in eachmatch(r"[+-]?\h*[i]?\h*[XYZI_]+", a)]
    Tableau(paulis)
end

macro T_str(a)
    quote _T_str($a) end
end

Base.getindex(tab::Tableau, i::Int) = PauliOperator(tab.phases[i], nqubits(tab), tab.xzs[:,i])
@inline function Base.getindex(tab::Tableau, r::Int, c::Int)
     # TODO this has code repetition with the Pauli getindex
     Tme = eltype(tab.xzs)
     x = (tab.xzs[_div(Tme,c-1)+1,r] & Tme(0x1)<<_mod(Tme,c-1))!=0x0
     z = (tab.xzs[end÷2+_div(Tme,c-1)+1,r] & Tme(0x1)<<_mod(Tme,c-1))!=0x0
     return (x,z)
end
Base.getindex(tab::Tableau, r) = Tableau(tab.phases[r], nqubits(tab), tab.xzs[:,r])
Base.getindex(tab::Tableau, r::Union{Colon,AbstractVector}, c::Union{Colon,AbstractVector}) = Tableau([s[c] for s in tab[r]])
Base.getindex(tab::Tableau, r::Union{Colon,AbstractVector}, c::Int) = tab[r,[c]]
Base.getindex(tab::Tableau, r::Int, c::Union{Colon,AbstractVector}) = tab[r][c]
Base.view(tab::Tableau, r) = Tableau(view(tab.phases, r), nqubits(tab), view(tab.xzs, :, r))

Base.iterate(tab::Tableau, state::Int=1) = state>length(tab) ? nothing : (tab[state], state+1)

function Base.setindex!(tab::Tableau, pauli::PauliOperator, i)
    tab.phases[i] = pauli.phase[]
    #tab.xzs[:,i] = pauli.xz # TODO why is this assigment causing allocations
    for j in 1:length(pauli.xz)
        tab.xzs[j,i] = pauli.xz[j]
    end
    tab
end

function Base.setindex!(tab::Tableau, t::Tableau, i)
    tab.phases[i] = t.phases
    tab.xzs[:,i] = t.xzs
    tab
end

function Base.setindex!(tab::Tableau{Tzv,Tm}, (x,z)::Tuple{Bool,Bool}, i, j) where {Tzv<:AbstractVector{UInt8}, Tme<:Unsigned, Tm<:AbstractMatrix{Tme}} # TODO this has code repetitions with the Pauli setindex
    if x
        tab.xzs[_div(Tme,j-1)+1,        i] |= Tme(0x1)<<_mod(Tme,j-1)
    else
        tab.xzs[_div(Tme,j-1)+1,        i] &= ~(Tme(0x1)<<_mod(Tme,j-1))
    end
    if z
        tab.xzs[end÷2+_div(Tme,j-1)+1, i] |= Tme(0x1)<<_mod(Tme,j-1)
    else
        tab.xzs[end÷2+_div(Tme,j-1)+1, i] &= ~(Tme(0x1)<<_mod(Tme,j-1))
    end
    tab
end

Base.firstindex(tab::Tableau) = 1

Base.lastindex(tab::Tableau) = length(tab.phases)
Base.lastindex(tab::Tableau, i) = size(tab)[i]

Base.eachindex(tab::Tableau) = Base.OneTo(lastindex(tab.phases))

Base.axes(tab::Tableau) = (Base.OneTo(lastindex(tab)), Base.OneTo(nqubits(tab)))
Base.axes(tab::Tableau,i) = axes(tab)[i]

Base.size(tab::Tableau) = (length(tab.phases),nqubits(tab))
Base.size(tab::Tableau,i) = size(tab)[i]

Base.length(tab::Tableau) = length(tab.phases)

function _show(io::IO, t::Tableau, limit=50, limit_vertical=20)
    padding = limit_vertical÷2-3
    n = size(t,1)
    range = 1:n
    if (limit_vertical < n && limit_vertical != -1)
        range = [1:padding; -1; (n-padding):n]
    end
    for i in range
        if (i == -1)
            print(io," ⋮\n")
            continue
        end
        _show(io, t[i], limit-7)
        i!=n && println(io)
    end
end

function Base.show(io::IO, t::Tableau)
    if get(io, :compact, false)
        r,q = size(t)
        print(io, "Tableaux $r×$q")
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, t, sz[2], sz[1])
    else
        _show(io, t, -1, -1)
    end
end

Base.:(==)(l::Tableau, r::Tableau) = r.nqubits==l.nqubits && r.phases==l.phases && r.xzs==l.xzs

Base.hash(t::Tableau, h::UInt) = hash(t.nqubits, hash(t.phases, hash(t.xzs, h)))

Base.copy(t::Tableau) = Tableau(copy(t.phases), t.nqubits, copy(t.xzs))

function Base.zero(::Type{Tableau{Tzv, Tm}}, r, q) where {Tzv,T<:Unsigned,Tm<:AbstractMatrix{T}}
    phases = zeros(UInt8,r)
    xzs = zeros(UInt, _nchunks(q,T), r)
    Tableau(phases, q, xzs)::Tableau{Vector{UInt8},Matrix{UInt}}
end
Base.zero(::Type{Tableau}, r, q) = zero(Tableau{Vector{UInt8},Matrix{UInt}}, r, q)
Base.zero(::Type{T}, q) where {T<:Tableau}= zero(T, q, q)
Base.zero(s::T) where {T<:Tableau} = zero(T, length(s), nqubits(s))

"""Zero-out a given row of a [`Tableau`](@ref)"""
@inline function zero!(t::Tableau,i)
    t.xzs[:,i] .= zero(eltype(t.xzs))
    t.phases[i] = 0x0
    t
end

##############################
# Stabilizer States
##############################

abstract type AbstractQCState end # This could include classical bits

abstract type AbstractStabilizer <: AbstractQCState end # This includes only qubits in stabilizer states

"""
Stabilizer, i.e. a list of commuting multi-qubit Hermitian Pauli operators.

Instances can be created with the `S` custom string macro or
as direct sum of other stabilizers.

!!! tip "Stabilizers and Destabilizers"
    In many cases you probably would prefer to use the [`MixedDestabilizer`](@ref)
    data structure, as it caries a lot of useful additional information, like tracking
    rank and destabilizer operators. `Stabilizer` has mostly a pedagogical value, and it
    is also used for slightly faster simulation of a particular subset of Clifford
    operations.

```jldoctest stabilizer
julia> s = S"XXX
             ZZI
             IZZ"
+ XXX
+ ZZ_
+ _ZZ

julia> s⊗s
+ XXX___
+ ZZ____
+ _ZZ___
+ ___XXX
+ ___ZZ_
+ ____ZZ
```

It has an indexing API, looking like a list of `PauliOperator`s.

```jldoctest stabilizer
julia> s[2]
+ ZZ_
```

Pauli operators can act directly on the a stabilizer.

```jldoctest stabilizer
julia> P"YYY" * s
- XXX
+ ZZ_
+ _ZZ
```

There are a number of ways to create a Stabilizer, including:

- generate Stabilizers from a list of Pauli operators

```jldoctest stabilizer
julia> Stabilizer([P"XX", P"ZZ"])
+ XX
+ ZZ
```

- generate Stabilizers from boolean matrices

```jldoctest stabilizer
julia> a = [true true; false false]; b = [false true; true true];

julia> Stabilizer(a, b)
+ XY
+ ZZ

julia> Stabilizer([0x0, 0x2], a, b)
+ XY
- ZZ
```

- initialize an empty Stabilizer and fill it through indexing

```jldoctest stabilizer
julia> s = zero(Stabilizer, 2)
+ __
+ __

julia> s[1,1] = (true, false); s
+ X_
+ __
```

There are no automatic checks for correctness (i.e. independence of all rows,
commutativity of all rows, hermiticity of all rows). The rank (number of rows)
is permitted to be less than the number of qubits (number of columns):
canonilization, projection, etc. continue working in that case. To great extent
this library uses the `Stabilizer` data structure simply as a tableau. This
might be properly abstracted away in future versions.

See also: [`PauliOperator`](@ref), [`canonicalize!`](@ref)
"""
struct Stabilizer{T<:Tableau} <: AbstractStabilizer
    tab::T
end

Stabilizer(phases::Tzv, nqubits::Int, xzs::Tm) where {Tzv<:AbstractVector{UInt8}, Tm<:AbstractMatrix{<:Unsigned}} = Stabilizer(Tableau(phases, nqubits, xzs))
Stabilizer(paulis::AbstractVector{PauliOperator{Tz,Tv}}) where {Tz,Tv} = Stabilizer(Tableau(paulis))
Stabilizer(phases::AbstractVector{UInt8}, xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Stabilizer(Tableau(phases, xs, zs))
Stabilizer(phases::AbstractVector{UInt8}, xzs::AbstractMatrix{Bool}) = Stabilizer(Tableau(phases, xzs))
Stabilizer(xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Stabilizer(Tableau(xs,zs))
Stabilizer(xzs::AbstractMatrix{Bool}) = Stabilizer(Tableau(xzs))
Stabilizer(s::Stabilizer) = s
macro S_str(a)
    quote Stabilizer(_T_str($a)) end
end
Base.getindex(stab::Stabilizer, i::Int) = tab(stab)[i]
Base.getindex(stab::Stabilizer, i) = Stabilizer(tab(stab)[i])
@inline Base.getindex(stab::Stabilizer, r::Int, c) = tab(stab)[r,c]
Base.getindex(stab::Stabilizer, r, c) = Stabilizer(tab(stab)[r,c])
Base.view(stab::Stabilizer, r) = Stabilizer(view(tab(stab),r))
Base.iterate(stab::Stabilizer, state=1) = iterate(tab(stab), state)
Base.setindex!(stab::Stabilizer, pauli::PauliOperator, i) = setindex!(tab(stab), pauli, i)
Base.setindex!(stab::Stabilizer, s::Stabilizer, i) = setindex!(tab(stab), tab(s), i)
Base.setindex!(stab::Stabilizer, (x,z)::Tuple{Bool,Bool}, i, j) = setindex!(tab(stab), (x,z), i, j)
Base.firstindex(stab::Stabilizer) = firstindex(tab(stab))
Base.lastindex(stab::Stabilizer) = lastindex(tab(stab))
Base.lastindex(stab::Stabilizer, i) = lastindex(tab(stab),i)
Base.eachindex(stab::Stabilizer) = eachindex(tab(stab))
Base.axes(stab::Stabilizer) = axes(tab(stab))
Base.axes(stab::Stabilizer, i) = axes(tab(stab), i)
Base.size(stab::Stabilizer) = size(tab(stab))
Base.size(stab::Stabilizer,i) = size(tab(stab),i)
Base.length(stab::Stabilizer) = length(tab(stab))

Base.:(==)(l::Stabilizer, r::Stabilizer) = tab(l) == tab(r)
Base.hash(s::Stabilizer, h::UInt) = hash(tab(s), h)
Base.copy(s::Stabilizer) = Stabilizer(copy(tab(s)))
Base.zero(::Type{S}, q) where {S<:Stabilizer} = zero(S, q, q)
Base.zero(::Type{Stabilizer{T}}, r, q) where {T<:Tableau} = Stabilizer(zero(T, r, q))
Base.zero(::Type{Stabilizer}, r, q) = zero(Stabilizer{Tableau{Vector{UInt8},Matrix{UInt}}}, r, q)
Base.zero(s::S) where {S<:Stabilizer} = zero(S, length(s), nqubits(s))
@inline zero!(s::Stabilizer,i) = zero!(tab(s),i)

function Base.show(io::IO, s::Stabilizer)
    if get(io, :compact, false)
        r,q = size(s)
        print(io, "Stabilizer $r×$q")
    else
        show(io, tab(s))
    end
end

##############################
# Helpers for subclasses of AbstractStabilizer that use Stabilizer as a tableau internally.
##############################

function Base.:(==)(l::T, r::S; phases=true) where {T<:AbstractStabilizer, S<:AbstractStabilizer}
    if phases
        return T==S && tab(l)==tab(r)
    else
        return T==S && tab(l).xzs==tab(r).xzs
    end
end

Base.hash(s::T, h::UInt) where {T<:AbstractStabilizer} = hash(T, hash(tab(s), h))

"""Extract the underlying tableau structure.

```jldoctest
julia> s = S"X"
+ X

julia> tab(s)
+ X

julia> tab(Destabilizer(s))
+ Z
+ X

julia> tab(MixedDestabilizer(s))
+ Z
+ X

julia> tab(tHadamard)
+ Z
+ X

julia> typeof(tab(tHadamard))
QuantumClifford.Tableau{Vector{UInt8}, Matrix{UInt64}}
```

See also: [`stabilizerview`](@ref), [`destabilizerview`](@ref), [`logicalxview`](@ref), [`logicalzview`](@ref)
"""
tab(s::Stabilizer{T}) where {T} = s.tab
tab(s::AbstractStabilizer) = s.tab

##############################
# Destabilizer formalism
##############################

"""
A tableau representation of a pure stabilizer state. The tableau tracks the
destabilizers as well, for efficient projections. On initialization there are
no checks that the provided state is indeed pure. This enables the use of this
data structure for mixed stabilizer state, but a better choice would be to use
[`MixedDestabilizer`](@ref).
""" # TODO clean up and document constructor
struct Destabilizer{T<:Tableau} <: AbstractStabilizer
    tab::T
end

function Destabilizer(s::Stabilizer)
    row, col = size(s)
    row>col && error(DomainError("The input stabilizer has more rows than columns, making it inconsistent or overdetermined."))
    mixed_destab = MixedDestabilizer(s)
    t = vcat(tab(destabilizerview(mixed_destab)),tab(stabilizerview(mixed_destab)))
    Destabilizer(t)
end

function Base.show(io::IO, d::Destabilizer)
    if get(io, :compact, false)
        r,q = size(stabilizerview(d))
        print(io, "Destablizer $r×$q")
    elseif get(io, :limit, false)
        h,w = displaysize(io)
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(min(w-9,size(d.tab,2)-4),0))
        _show(io, destabilizerview(d).tab, w, h÷2)
        println(io, "\n𝒮𝓉𝒶𝒷" * "━"^max(min(w-7,size(d.tab,2)-2),0))
        _show(io, stabilizerview(d).tab, w, h÷2)
    else
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-4,0))
        _show(io, destabilizerview(d).tab, -1, -1)
        println(io, "\n𝒮𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-2,0))
        _show(io, stabilizerview(d).tab, -1, -1)
    end
end

Base.length(d::Destabilizer) = length(tab(d))÷2

Base.copy(d::Destabilizer) = Destabilizer(copy(tab(d)))

##############################
# Mixed Stabilizer states
##############################

"""
A slight improvement of the [`Stabilizer`](@ref) data structure that enables
more naturally and completely the treatment of mixed states, in particular when
the [`project!`](@ref) function is used.
"""
mutable struct MixedStabilizer{T<:Tableau} <: AbstractStabilizer
    tab::T # TODO assert size on construction
    rank::Int
end

function MixedStabilizer(s::Stabilizer{T}) where {T}
    s, xr, zr = canonicalize!(s,ranks=true)
    spadded = zero(Stabilizer, nqubits(s))
    spadded[1:zr] = s[1:zr]
    MixedStabilizer(spadded,zr)
end

MixedStabilizer(s::Stabilizer,rank::Int) = MixedStabilizer(tab(s),rank)

function Base.show(io::IO, s::MixedStabilizer)
    if get(io, :compact, false)
        r,q = size(s)
        print(io, "MixedStabilizer $r×$q")
    else
        show(io, stabilizerview(s))
    end
end

Base.length(d::MixedStabilizer) = length(d.tab)

Base.copy(ms::MixedStabilizer) = MixedStabilizer(copy(ms.tab),ms.rank)

##############################
# Mixed Destabilizer states
##############################

"""
A tableau representation for mixed stabilizer states that keeps track of the
destabilizers in order to provide efficient projection operations.

The rank `r` of the `n`-qubit tableau is tracked, either so that it can be
used to represent a mixed stabilizer state, or so that it can be used to
represent an `n-r` logical-qubit code over `n` physical qubits. The "logical"
operators are tracked as well.

When the constructor is called on an incomplete [`Stabilizer`](@ref) it
automatically calculates the destabilizers and logical operators, following
chapter 4 of [gottesman1997stabilizer](@cite).

See also: [`stabilizerview`](@ref), [`destabilizerview`](@ref), [`logicalxview`](@ref), [`logicalzview`](@ref)
"""
mutable struct MixedDestabilizer{T<:Tableau} <: AbstractStabilizer
    tab::T # TODO assert size on construction
    rank::Int
end

# Added a lot of type assertions to help Julia infer types
function MixedDestabilizer(stab::Stabilizer{T}; undoperm=true) where {T}
    rows,n = size(stab)
    stab, r, s, permx, permz, ops = canonicalize_gott!(copy(stab))
    t = zero(T, n*2, n)
    vstab = @view tab(stab)[1:r+s] # this view is necessary for cases of tableaux with redundant rows
    t[n+1:n+r+s] = vstab # The Stabilizer part of the tableau
    for i in 1:r # The Destabilizer part
        t[i,i] = (false,true)
    end
    for i in r+1:r+s
        t[i,i] = (true,false)
    end
    if r+s!=n
        H = stab_to_gf2(vstab)            # n-k × 2n
        k = n - r - s
        E = H[r+1:end,end÷2+r+s+1:end]    # n-k-r × n-r-s
        C1 = H[1:r,end÷2+r+1:end÷2+r+s]   #     r ×     s
        C2 = H[1:r,end÷2+r+s+1:end]       #     r × n-r-s
        U2 = E'                           # n-r-s × n-k-r
        V1 = (E' * C1' + C2').%2 .!= 0x0  # n-r-s ×     r
        #X = hcat(zeros(Bool,k,r),U2,I,V1,zeros(Bool,k,s+k))
        X = zeros(Bool, k, 2n)
        X[:,r+1:n-k] .= U2
        X[:,n+1:n+r] .= V1
        @inbounds @simd for i in 1:k X[i,n-k+i] = true end
        sX = Tableau(X)
        t[r+s+1:n] = sX # One of the logical sets in the tableau
        A2 = H[1:r,r+s+1:end÷2]           #     r × n-r-s
        #Z = hcat(zeros(Bool,k,n),A2',zeros(Bool,k,s),I)
        Z = zeros(Bool, k, 2n)
        Z[:,n+1:n+r] .= A2'
        @inbounds @simd for i in 1:k Z[i,2n-k+i] = true end
        sZ = Tableau(Z)
        t[n+r+s+1:end] = sZ # The other logical set in the tableau
    end
    if undoperm
        # print("s, r = ", s, r)
        for op in reverse(ops)
            # print("\n\nPerforming: ", op, '\n')
            if op[1] == 1
                # show(IOContext(stdout::IO, :limit => true), t)   
                # print('\n')
                # print("=>\n")

                rowswap!(t, op[2], op[3]; phases=Val(true))
                rowswap!(t, op[2] + r + s, op[3] + r + s; phases=Val(true))
                # show(IOContext(stdout::IO, :limit => true), t)   
                # print('\n')

            elseif op[1] == 2
                # if the mul_left ops do not need be reversed (to keep linear independence),
                # i will take the old approach of just modifying permx, permz

                # show(IOContext(stdout::IO, :limit => true), t)   
                # print('\n')
                # print("=>\n")

                # mul_left!(t, op[2], op[3]; phases=Val(true))
                # mul_left!(t, op[2] + r + s, op[3] + r + s; phases=Val(true))

                # mul_left!(t, op[2], op[3]; phases=Val(true))
                # mul_left!(t, op[2] + r + s, op[3] + r + s; phases=Val(true))

                # mul_left!(t, op[2], op[3]; phases=Val(true))
                # mul_left!(t, op[2] + r + s, op[3] + r + s; phases=Val(true))

                # show(IOContext(stdout::IO, :limit => true), t)   
                # print('\n')

            elseif op[1] == 3
                t = t[:,invperm(permx)]::T
            elseif op[1] == 4
                t = t[:,invperm(permz)]::T
            end
        end
    end
    MixedDestabilizer(t, r+s)::MixedDestabilizer{T}
end

function MixedDestabilizer(d::Destabilizer, r::Int)
    l,n = size(d.tab)
    if l==2n
        MixedDestabilizer(tab(d), r)
    else
        throw(DomainError("Only full-rank `Destabilizer` can be converted to `MixedDestabilizer` with specific rank. Try not specifying `r`."))
    end
end
function MixedDestabilizer(d::Destabilizer)
    l,n = size(d.tab)
    if l==2n
        MixedDestabilizer(d, nqubits(d))
    else
        MixedDestabilizer(stabilizerview(d))
    end
end

function MixedDestabilizer(d::MixedStabilizer) MixedDestabilizer(stabilizerview(d)) end

function Base.show(io::IO, d::MixedDestabilizer)
    r = rank(d)
    q = nqubits(d)
    if get(io, :compact, false)
        print(io, "MixedDestablizer $r×$q")
    elseif get(io, :limit, false)
        h,w = displaysize(io)
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(min(w-9,size(d.tab,2)-4),0))
        _show(io, destabilizerview(d).tab, w, h÷4)
        if r != q
            println(io)
            println(io, "𝒳ₗ" * "━"^max(min(w-5,size(d.tab,2)),0))
            _show(io, logicalxview(d).tab, w, h÷4)
        end
        println(io)
        println(io, "𝒮𝓉𝒶𝒷" * "━"^max(min(w-7,size(d.tab,2)-2),0))
        _show(io, stabilizerview(d).tab, w, h÷4)
        if r != q
            println(io)
            println(io, "𝒵ₗ" * "━"^max(min(w-5,size(d.tab,2)),0))
            _show(io, logicalzview(d).tab, w, h÷4)
        end
    else
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-4,0))
        _show(io, destabilizerview(d).tab, -1, -1)
        if r != q
            println(io)
            println(io, "𝒳ₗ" * "━"^max(size(d.tab,2),0))
            _show(io, logicalxview(d).tab, -1, -1)
        end
        println(io)
        println(io, "𝒮𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-2,0))
        _show(io, stabilizerview(d).tab, -1, -1)
        if r != q
            println(io)
            println(io, "𝒵ₗ" * "━"^max(size(d.tab,2)),0)
            _show(io, logicalzview(d).tab, -1, -1)
        end
    end
end

Base.length(d::MixedDestabilizer) = length(d.tab)÷2

Base.copy(d::MixedDestabilizer) = MixedDestabilizer(copy(d.tab),d.rank)

##############################
# Subtableau views
##############################

"""A view of the subtableau corresponding to the stabilizer. See also [`tab`](@ref), [`destabilizerview`](@ref), [`logicalxview`](@ref), [`logicalzview`](@ref)"""
@inline stabilizerview(s::Stabilizer) = s
@inline stabilizerview(s::Destabilizer) = Stabilizer(@view tab(s)[end÷2+1:end])
@inline stabilizerview(s::MixedStabilizer) = Stabilizer(@view tab(s)[1:s.rank])
@inline stabilizerview(s::MixedDestabilizer) = Stabilizer(@view tab(s)[end÷2+1:end÷2+s.rank])

"""A view of the subtableau corresponding to the destabilizer. See also [`tab`](@ref), [`stabilizerview`](@ref), [`logicalxview`](@ref), [`logicalzview`](@ref)"""
@inline destabilizerview(s::Destabilizer) = Stabilizer(@view tab(s)[1:end÷2])
@inline destabilizerview(s::MixedDestabilizer) = Stabilizer(@view tab(s)[1:s.rank])

"""A view of the subtableau corresponding to the logical X operators. See also [`tab`](@ref), [`stabilizerview`](@ref), [`destabilizerview`](@ref), [`logicalzview`](@ref)"""
@inline logicalxview(s::MixedDestabilizer) = Stabilizer(@view tab(s)[s.rank+1:end÷2])
"""A view of the subtableau corresponding to the logical Z operators. See also [`tab`](@ref), [`stabilizerview`](@ref), [`destabilizerview`](@ref), [`logicalxview`](@ref)"""
@inline logicalzview(s::MixedDestabilizer) = Stabilizer(@view tab(s)[end÷2+s.rank+1:end])

"""The number of qubits of a given state."""
@inline nqubits(s::AbstractStabilizer) = nqubits(tab(s))
@inline nqubits(t::Tableau) = t.nqubits

"""The phases of a given tableau."""
@inline phases(t::Tableau) = t.phases
@inline phases(s::Stabilizer) = phases(tab(s))

##############################
# Pauli Operator Helpers
##############################

"""
Get the phase of the product of two Pauli operators.

Phase is encoded as F(4) in the low qubits of an UInt8.

```jldoctest
julia> P"ZZZ"*P"XXX"
-iYYY

julia> prodphase(P"ZZZ", P"XXX")
0x03

julia> prodphase(P"XXX", P"ZZZ")
0x01
```
"""
@inline function prodphase(l::AbstractVector{T}, r::AbstractVector{T})::UInt8 where T<:Unsigned
    res = 0
    len = length(l)÷2
    @inbounds @simd for i in 1:len
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        res += count_ones((~z2 & x2 & ~x1 & z1) | ( z2 & ~x2 & x1 & z1) | (z2 &  x2 & x1 & ~z1))
        res -= count_ones(( z2 & x2 & ~x1 & z1) | (~z2 &  x2 & x1 & z1) | (z2 & ~x2 & x1 & ~z1))
    end
    unsigned(res)&0x3
end

"""The quantumlib/Stim implementation, which performs the prodphase and mul_left! together. Used for unit tests."""
function _stim_prodphase(l::AbstractVector{T}, r::AbstractVector{T}) where T<: Unsigned
    cnt1 = zero(T)
    cnt2 = zero(T)
    len = length(l)÷2
    @inbounds @simd for i in 1:len
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        newx1 = x1 ⊻ x2 # Here l or r would be updated in an actual multiplication
        newz1 = z1 ⊻ z2 # Here l or r would be updated in an actual multiplication
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        cnt2 ⊻= (cnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        cnt1 ⊻= anti_comm
    end
    s = count_ones(cnt1)
    s ⊻= count_ones(cnt2) << 1
    s & 0x3
end

@inline function prodphase(l::PauliOperator, r::PauliOperator)::UInt8
    (l.phase[]+r.phase[]+prodphase(l.xz,r.xz))&0x3
end

@inline function prodphase(l::PauliOperator, r::Tableau, i)::UInt8
    (l.phase[]+r.phases[i]+prodphase(l.xz, (@view r.xzs[:,i])))&0x3
end

@inline function prodphase(l::Tableau, r::PauliOperator, i)::UInt8
    (l.phases[i]+r.phase[]+prodphase((@view l.xzs[:,i]), r.xz))&0x3
end

@inline function prodphase(l::Tableau, r::Tableau, i, j)::UInt8
    (l.phases[i]+r.phases[j]+prodphase((@view l.xzs[:,i]), (@view r.xzs[:,j])))&0x3
end

@inline prodphase(l::PauliOperator,r::Stabilizer,i) = prodphase(l,tab(r),i)
@inline prodphase(l::Stabilizer,r::PauliOperator,i) = prodphase(tab(l),r,i)
@inline prodphase(l::Stabilizer,r::Stabilizer,i,j) = prodphase(tab(l),tab(r),i,j)

"""
Check whether two operators commute.

`0x0` if they commute, `0x1` if they anticommute.

```jldoctest
julia> P"XX"*P"ZZ", P"ZZ"*P"XX"
(- YY, - YY)

julia> comm(P"ZZ", P"XX")
0x00

julia> comm(P"IZ", P"XX")
0x01
```
"""
@inline function comm(l::AbstractVector{T}, r::AbstractVector{T})::UInt8 where T<:Unsigned
    res = T(0)
    len = length(l)÷2
    @inbounds @simd for i in 1:len
        res ⊻= (l[i+len] & r[i]) ⊻ (l[i] & r[i+len])
    end
    count_ones(res)&0x1
end

@inline function comm(l::PauliOperator, r::PauliOperator)::UInt8
    comm(l.xz,r.xz)
end

@inline function comm(l::PauliOperator, r::Tableau, i::Int)::UInt8
    comm(l.xz,(@view r.xzs[:,i]))
end

@inline function comm(l::Tableau, r::PauliOperator, i::Int)::UInt8
    comm(r, l, i)
end

@inline function comm(s::Tableau, l::Int, r::Int)::UInt8
    comm((@view s.xzs[:,l]),(@view s.xzs[:,r]))
end

function comm(l::PauliOperator, r::Tableau)::Vector{UInt8}
    [comm(l,r,i) for i in 1:size(r,1)]
end

comm(l::Tableau, r::PauliOperator) = comm(r, l)

@inline comm(l::PauliOperator, r::Stabilizer, i::Int) = comm(l, tab(r), i)
@inline comm(l::Stabilizer, r::PauliOperator, i::Int) = comm(tab(l), r, i)
@inline comm(l::PauliOperator, r::Stabilizer) = comm(l, tab(r))
@inline comm(l::Stabilizer, r::PauliOperator) = comm(tab(l), r)
@inline comm(s::Stabilizer, l::Int, r::Int) = comm(tab(s), l, r)

Base.:(*)(l::PauliOperator, r::PauliOperator) = mul_left!(copy(r),l)

"""Nonvectorized version of `mul_left!` used for unit tests."""
function _mul_left_nonvec!(r::AbstractVector{T}, l::AbstractVector{T}; phases::Bool=true) where T<:Unsigned
    if !phases
        r .⊻= l
        return zero(T)
    end
    cnt1 = zero(T)
    cnt2 = zero(T)
    len = length(l)÷2
    @inbounds @simd for i in 1:len
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        r[i] = newx1 = x1 ⊻ x2
        r[i+len] = newz1 = z1 ⊻ z2
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        cnt2 ⊻= (cnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        cnt1 ⊻= anti_comm
    end
    s = count_ones(cnt1)
    s ⊻= count_ones(cnt2) << 1
    s
end

function mul_left!(r::AbstractVector{T}, l::AbstractVector{T}; phases::Val{B}=Val(true))::UInt8 where {T<:Unsigned, B}
    if !B
        r .⊻= l
        return UInt8(0x0)
    end
    len = length(l)÷2
    veclen = Int(pick_vector_width(T)) # TODO remove the Int cast
    rcnt1 = 0
    rcnt2 = 0
    packs = len÷veclen
    VT = SIMD.Vec{veclen,T}
    if packs>0
        cnt1 = zero(VT)
        cnt2 = zero(VT)
        lane = SIMD.VecRange{veclen}(0)
        @inbounds for i in 1:veclen:(len-veclen+1)
            # JET-XXX The ::VT should not be necessary, but they help with inference
            x1::VT, x2::VT, z1::VT, z2::VT = l[i+lane], r[i+lane], l[i+len+lane], r[i+len+lane]
            r[i+lane] = newx1 = x1 ⊻ x2
            r[i+len+lane] = newz1 = z1 ⊻ z2
            x1z2 = x1 & z2
            anti_comm = (x2 & z1) ⊻ x1z2
            cnt2 ⊻= (newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
            cnt1 ⊻= anti_comm
        end
        for i in 1:length(cnt1)
            rcnt1 += count_ones(cnt1[i])
        end
        for i in 1:length(cnt2)
            rcnt2 += count_ones(cnt2[i])
        end
    end
    scnt1 = zero(T)
    scnt2 = zero(T)
    @inbounds for i in (packs*veclen+1):len # same code for the pieces that do not fit in a vector at the end (TODO padding would simplify the code)
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        r[i] = newx1 = x1 ⊻ x2
        r[i+len] = newz1 = z1 ⊻ z2
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        scnt2 ⊻= (scnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        scnt1 ⊻= anti_comm
    end
    rcnt1 += count_ones(scnt1)
    rcnt2 += count_ones(scnt2)
    UInt8((rcnt1 ⊻ (rcnt2<<1))&0x3)
end

@inline function mul_left!(r::PauliOperator, l::PauliOperator; phases::Val{B}=Val(true)) where B
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The two Pauli operators should have the same length!")) # TODO skip this when @inbounds is set
    s = mul_left!(r.xz, l.xz, phases=phases)
    B && (r.phase[] = (s+r.phase[]+l.phase[])&0x3)
    r
end

@inline function mul_left!(r::PauliOperator, l::Tableau, i::Int; phases::Val{B}=Val(true)) where B
    s = mul_left!(r.xz, (@view l.xzs[:,i]), phases=phases)
    B && (r.phase[] = (s+r.phase[]+l.phases[i])&0x3)
    r
end

@inline mul_left!(r::PauliOperator, l::Stabilizer, i::Int; phases::Val{B}=Val(true)) where B = mul_left!(r, tab(l), i; phases=phases)

(⊗)(l::PauliOperator, r::PauliOperator) = PauliOperator((l.phase[]+r.phase[])&0x3, vcat(xbit(l),xbit(r)), vcat(zbit(l),zbit(r)))

function Base.:(*)(l::Number, r::PauliOperator)
    p = copy(r)
    if l==1
        nothing
    elseif l==1im
        p.phase[] = (p.phase[] + 1)&0x3
    elseif l==-1
        p.phase[] = (p.phase[] + 2)&0x3
    elseif l==-1im
        p.phase[] = (p.phase[] + 3)&0x3
    else
        throw(DomainError(l,"Only {±1,±i} are permitted as phases."))
    end
    p
end

Base.:(+)(p::PauliOperator) = p

function Base.:(-)(p::PauliOperator)
    p = copy(p)
    p.phase[] = (p.phase[]+2)&0x3
    p
end

const I = P"I"
const Z = P"Z"
const X = P"X"
const Y = P"Y"

##############################
# Stabilizer helpers
##############################

@inline function rowswap!(s::Tableau, i, j; phases::Val{B}=Val(true)) where B # Written only so we can avoid copying in `canonicalize!`
    (i == j) && return
    B && begin s.phases[i], s.phases[j] = s.phases[j], s.phases[i] end
    @inbounds @simd for k in 1:size(s.xzs,1)
        s.xzs[k,i], s.xzs[k,j] = s.xzs[k,j], s.xzs[k,i]
    end
end

@inline function rowswap!(s::Stabilizer, i, j; phases::Val{B}=Val(true)) where B
    rowswap!(tab(s), i, j; phases)
end

@inline function rowswap!(s::Destabilizer, i, j; phases::Val{B}=Val(true)) where B
    t = tab(s)
    rowswap!(t, i, j; phases)
    n = size(t,1)÷2
    rowswap!(t, i+n, j+n; phases)
end

@inline function rowswap!(s::MixedStabilizer, i, j; phases::Val{B}=Val(true)) where B
    rowswap!(tab(s), i, j; phases)
end

@inline function rowswap!(s::MixedDestabilizer, i, j; phases::Val{B}=Val(true)) where B
    t = tab(s)
    rowswap!(t, i, j; phases)
    n = nqubits(s)
    rowswap!(t, i+n, j+n; phases)
end

@inline function mul_left!(s::Tableau, m, i; phases::Val{B}=Val(true)) where B
    extra_phase = mul_left!((@view s.xzs[:,m]), (@view s.xzs[:,i]); phases=phases)
    B && (s.phases[m] = (extra_phase+s.phases[m]+s.phases[i])&0x3)
    s
end

@inline function mul_left!(s::Stabilizer, m, i; phases::Val{B}=Val(true)) where B
    mul_left!(tab(s), m, i; phases)
end

@inline function mul_left!(s::Destabilizer, i, j; phases::Val{B}=Val(true)) where B
    t = tab(s)
    mul_left!(t, j, i; phases=Val(false)) # Indices are flipped to preserve commutation constraints
    n = size(t,1)÷2
    mul_left!(t, i+n, j+n; phases=phases)
end

@inline function mul_left!(s::MixedStabilizer, i, j; phases::Val{B}=Val(true)) where B
    mul_left!(tab(s), i, j; phases=phases)
end

@inline function mul_left!(s::MixedDestabilizer, i, j; phases::Val{B}=Val(true)) where B
    t = tab(s)
    mul_left!(t, j, i; phases=Val(false)) # Indices are flipped to preserve commutation constraints
    n = nqubits(t)
    mul_left!(t, i+n, j+n; phases=phases)
end

# TODO is this used anywhere?
"""Swap two columns of a stabilizer in place."""
@inline function colswap!(s::Tableau, i, j)
    Tme = eltype(s.xzs)
    lowbit = Tme(1)
    ibig = _div(Tme,i-1)+1
    ismall = _mod(Tme,i-1)
    ismallm = lowbit<<(ismall)
    jbig = _div(Tme,j-1)+1
    jsmall = _mod(Tme,j-1)
    jsmallm = lowbit<<(jsmall)
    for off in [0,size(s.xzs,2)÷2]
        ibig += off
        jbig += off
        @inbounds for k in 1:size(s.xzs,2)
            ival = s.xzs[ibig,k] & ismallm
            jval = s.xzs[jbig,k] & jsmallm
            s.xzs[ibig,k] &= ~ismallm
            s.xzs[jbig,k] &= ~jsmallm
            s.xzs[ibig,k] |= jval>>>(jsmall-ismall)
            s.xzs[jbig,k] |= ival>>>(ismall-jsmall)
        end
    end
end

@inline _logsizeof(::UInt128) = 7
@inline _logsizeof(::UInt64 ) = 6
@inline _logsizeof(::UInt32 ) = 5
@inline _logsizeof(::UInt16 ) = 4
@inline _logsizeof(::UInt8  ) = 3
@inline _logsizeof(::Type{UInt128}) = 7
@inline _logsizeof(::Type{UInt64 }) = 6
@inline _logsizeof(::Type{UInt32 }) = 5
@inline _logsizeof(::Type{UInt16 }) = 4
@inline _logsizeof(::Type{UInt8  }) = 3
@inline _mask(::T) where T<:Unsigned = sizeof(T)*8-1
@inline _mask(arg::T) where T<:Type = sizeof(arg)*8-1
@inline _div(T,l) = l >> _logsizeof(T)
@inline _mod(T,l) = l & _mask(T)

function unsafe_bitfindnext_(chunks::AbstractVector{T}, start::Int) where T<:Unsigned
    chunk_start = _div(T,start-1)+1
    within_chunk_start = _mod(T,start-1)
    mask = ~T(0) << within_chunk_start

    @inbounds begin
        if chunks[chunk_start] & mask != 0
            return (chunk_start-1) << 6 + trailing_zeros(chunks[chunk_start] & mask) + 1
        end

        for i = chunk_start+1:length(chunks)
            if chunks[i] != 0
                return (i-1) << 6 + trailing_zeros(chunks[i]) + 1
            end
        end
    end
    return nothing
end

"""Permute the qubits (i.e., columns) of the tableau in place."""
function Base.permute!(s::Tableau, perm::AbstractVector)
    for r in 1:size(s,1)
        s[r] = s[r][perm] # TODO make a local temporary buffer row instead of constantly allocating new rows
    end
    s
end

"""Permute the qubits (i.e., columns) of the state in place."""
function Base.permute!(s::AbstractStabilizer, perm::AbstractVector)
    permute!(tab(s), perm)
    s
end

function check_allrowscommute(stabilizer::Tableau)
    for i in eachindex(stabilizer)
        for j in eachindex(stabilizer)
            i==j && continue
            comm(stabilizer[i],stabilizer[j])==0x0 || return false
        end
    end
    return true
end
check_allrowscommute(stabilizer::Stabilizer)=check_allrowscommute(tab(stabilizer))

function Base.vcat(tabs::Tableau...)
    Tableau(
        vcat((s.phases for s in tabs)...),
        tabs[1].nqubits,
        hcat((s.xzs for s in tabs)...))
end

Base.vcat(stabs::Stabilizer...) = Stabilizer(vcat((tab(s) for s in stabs)...))

##############################
# Unitary Clifford Operations
##############################

function Base.:(*)(p::AbstractCliffordOperator, s::AbstractStabilizer; phases::Bool=true)
    s = copy(s)
    @valbooldispatch _apply!(s,p; phases=Val(phases)) phases
end
function apply!(stab::AbstractStabilizer, op::AbstractCliffordOperator; phases::Bool=true)
    @valbooldispatch _apply!(stab,op; phases=Val(phases)) phases
end
function apply!(stab::AbstractStabilizer, op::AbstractCliffordOperator, indices; phases::Bool=true)
    @valbooldispatch _apply!(stab,op,indices; phases=Val(phases)) phases
end

# TODO no need to track phases outside of stabview
function _apply!(stab::AbstractStabilizer, p::PauliOperator; phases::Val{B}=Val(true)) where B
    nqubits(stab)==nqubits(p) || throw(DimensionMismatch("The tableau and the Pauli operator need to act on the same number of qubits. Consider specifying an array of indices as a third argument to the `apply!` function to avoid this error."))
    s = tab(stab)
    B || return stab
    for i in eachindex(s)
        s.phases[i] = (s.phases[i]+comm(p,s,i)<<1+p.phase[]<<1)&0x3
    end
    stab
end

function _apply!(stab::AbstractStabilizer, p::PauliOperator, indices; phases::Val{B}=Val(true)) where B
    s = tab(stab)
    B || return stab
    newp = zero(typeof(p),nqubits(s)) # TODO this is an unnecessarily slow operation for something that can be sparse
    for (ii,i) in enumerate(indices)
        newp[i] = p[ii]
    end
    newp.phase .= p.phase
    for i in eachindex(s)
        s.phases[i] = (s.phases[i]+comm(newp,s,i)<<1+newp.phase[]<<1)&0x3
    end
    stab
end

##############################
# Helpers for binary codes
##############################

"""The F(2,2) matrix of a given tableau, represented as the concatenation of two binary matrices, one for X and one for Z."""
function stab_to_gf2(s::Tableau)::Matrix{Bool}
    r, n = size(s)
    H = zeros(Bool,r,2n)
    for iᵣ in 1:r
        @inbounds @simd for iₙ in 1:n
            H[iᵣ,iₙ], H[iᵣ,iₙ+n] = s[iᵣ,iₙ]
        end
    end
    H
end
stab_to_gf2(s::Stabilizer) = stab_to_gf2(tab(s))

"""Gaussian elimination over the binary field."""
function gf2_gausselim!(H)
    rows, cols = size(H)
    j = 1
    for c in 1:cols
        j>rows && break
        goodrow = findfirst(H[j:end,c])
        goodrow===nothing && continue
        goodrow += j-1
        @inbounds for d in 1:cols
            H[goodrow,d], H[j,d] = H[j,d], H[goodrow,d]
        end
        for r in 1:rows
            r==j && continue
            if H[r,c]
                H[r,:] .⊻= H[j,:]
            end
        end
        j += 1
    end
    H
end

# TODO check whether Nemo is faster in doing this and other gf2 methods.
"""Check whether a binary matrix is invertible."""
function gf2_isinvertible(H) # TODO can be smarter and exit earlier. And should check for squareness.
    ut = gf2_gausselim!(copy(H))
    all((ut[i, i] for i in 1:size(ut,1)))::Bool
end

"""Invert a binary matrix."""
function gf2_invert(H)
    id = zero(H)
    s = size(H,1)
    for i in 1:s id[i,i]=true end
    M = hcat(H,id)
    gf2_gausselim!(M)
    M[:,s+1:end]
end

"""The permutation of columns which turns a binary matrix into standard form. It is assumed the matrix has already undergone Gaussian elimination."""
function gf2_H_standard_form_indices(H)
    rows, cols = size(H)
    goodindices = Int[]
    j = 1
    r = 1
    for r in 1:rows
        i = findfirst(H[r,:])
        (i ∈ goodindices)::Bool && continue
        push!(goodindices, i)
    end
    badindices = [r for r in 1:goodindices[end] if !(r ∈ goodindices)]
    return vcat(goodindices, badindices, goodindices[end]+1:cols)
end

"""For a given F(2,2) parity check matrix, return the generator matrix."""
function gf2_H_to_G(H)
    # XXX it assumes that H is upper triangular (Gauss elimination, canonicalized, etc)
    # XXX careful, this is the binary code matrix - for the F(2,2) code you need to swap the x and z parts
    rows, cols = size(H)
    sindx = gf2_H_standard_form_indices(H)
    H = H[:,sindx]
    P = H[:,rows+1:end]'
    I = falses(cols-rows,cols-rows)
    for i in 1:size(I,1)
        I[i,i]=true
    end
    G = hcat(P,I)
    G[:,invperm(sindx)]
end

##############################
# Error classes
##############################
struct BadDataStructure <: Exception
    message::String
    whichop::Symbol
    whichstructure::Symbol
end

##############################
# Consistency checks
##############################

"""
Check basic consistency requirements of a stabilizer. Used in tests.
"""
function stab_looks_good(s)
    c = tab(canonicalize!(copy(s)))
    nrows, ncols = size(c)
    all((c.phases .== 0x0) .| (c.phases .== 0x2)) || return false
    H = stab_to_gf2(c)
    good_indices = reduce(|,H,dims=(1,))
    good_indices = good_indices[1:end÷2] .| good_indices[end÷2+1:end]
    colsok = ncols>nrows || all(good_indices) # TODO, this can be stricter
    colsok || return false
    good_indices = reduce(|,H,dims=(2,))
    rowsok = nrows>ncols || all(good_indices) # TODO, this can be stricter
    rowsok || return false
    check_allrowscommute(c) || return false
    return true
end

"""
Check basic consistency requirements of a mixed stabilizer. Used in tests.
"""
function mixed_stab_looks_good(s)
    s = stabilizerview(s)
    c = tab(canonicalize!(copy(s)))
    all((c.phases .== 0x0) .| (c.phases .== 0x2)) || return false
    H = stab_to_gf2(c)
    # Unlike for stabilizers, we do not expect all columns to have non-identity ops
    non_trivial_rows = reduce(|,H,dims=(2,))
    all(non_trivial_rows) || return false
    check_allrowscommute(c) || return false
    return true
end

"""
Check basic consistency requirements of a destabilizer. Used in tests.
"""
function destab_looks_good(destabilizer)
    s = stabilizerview(destabilizer)
    d = destabilizerview(destabilizer)
    stab_looks_good(s) || return false
    stab_looks_good(d) || return false
    for i in eachindex(s)
        comm(s[i],d[i])==0x1 || return false
        for j in eachindex(s)
            j==i && continue
            comm(s[i],d[j])==0x0 || return false
        end
    end
    return true
end

"""
Check basic consistency requirements of a mixed destabilizer. Used in tests.
"""
function mixed_destab_looks_good(destabilizer)
    s = stabilizerview(destabilizer)
    d = destabilizerview(destabilizer)
    x = logicalxview(destabilizer)
    z = logicalzview(destabilizer)
    check_allrowscommute(s) || return false
    mixed_stab_looks_good(s) || return false
    mixed_stab_looks_good(d) || return false
    mixed_stab_looks_good(x) || return false
    mixed_stab_looks_good(z) || return false
    for i in eachindex(s)
        comm(s[i],d[i])==0x1 || return false
        for j in eachindex(s)
            j==i && continue
            comm(s[i],d[j])==0x0 || return false
        end
        for j in eachindex(x)
            comm(s[i],x[j])==0x0 || return false
            comm(s[i],z[j])==0x0 || return false
        end
    end
    for i in eachindex(x)
        for j in eachindex(x)
            comm(x[i],x[j])==0x0 || return false
            comm(z[i],z[j])==0x0 || return false
            if i==j
                comm(x[i],z[j])==0x1 || return false
            else
                comm(x[i],z[j])==0x0 || return false
            end
        end
    end
    return true
end

include("canonicalization.jl")
include("dense_cliffords.jl")
include("project_trace_reset.jl")
include("linalg.jl")
include("symbolic_cliffords.jl")
include("mctrajectory.jl")
include("misc_ops.jl")
include("classical_register.jl")
include("enumeration.jl")
include("randoms.jl")
include("useful_states.jl")
include("noise.jl")
include("affectedqubits.jl")
include("experimental/Experimental.jl")
include("graphs.jl")
include("entanglement.jl")
include("precompiles.jl")

end #module
