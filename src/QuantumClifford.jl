"""
A module for using the Stabilizer formalism and simulating Clifford circuits.
"""
module QuantumClifford

# TODO document phases=false

# TODO Significant performance improvements: many operations do not need phase=true if the Pauli operations commute

import LinearAlgebra
using DocStringExtensions
using Polyester
#using LoopVectorization
using HostCPUFeatures: pick_vector_width
import SIMD

export @P_str, PauliOperator, ⊗, I, X, Y, Z, permute,
    xbit, zbit, xview, zview,
    @S_str, Stabilizer, prodphase, comm, check_allrowscommute,
    Destabilizer, MixedStabilizer, MixedDestabilizer,
    nqubits, stabilizerview, destabilizerview, logicalxview, logicalzview,
    bitview, quantumstate,
    canonicalize!, canonicalize_rref!, canonicalize_gott!, colpermute!,
    logdot, expect,
    generate!, project!, reset_qubits!, traceout!,
    apply!,
    projectX!, projectY!, projectZ!,
    projectrand!, projectXrand!, projectYrand!, projectZrand!,
    tab, puttableau!,
    CliffordOperator, @C_str,
    tCNOT, tCPHASE, tSWAP, tHadamard, tPhase, tId1,
    AbstractSymbolicOperator, AbstractSingleQubitOperator, AbstractTwoQubitOperator,
    sHadamard, sPhase, sInvPhase, SingleQubitOperator, sId1, sX, sY, sZ,
    sCNOT, sCPHASE, sSWAP,
    SparseGate,
    sMX, sMY, sMZ, PauliMeasurement, Reset,
    Register,
    tensor, tensor_pow,
    stab_to_gf2, gf2_gausselim!, gf2_isinvertible, gf2_invert, gf2_H_to_G,
    single_z, single_x, single_y,
    enumerate_single_qubit_gates, random_clifford1,
    enumerate_cliffords, symplecticGS, clifford_cardinality, enumerate_phases,
    random_invertible_gf2,
    random_pauli, random_stabilizer, random_destabilizer, random_clifford,
    bell, ghz,
    BadDataStructure,
    graphstate, graphstate!, graph_gatesequence, graph_gate,
    canonicalize_clip!, bigram, entanglement_entropy


const BIG_INT_MINUS_ONE = Ref{BigInt}()
const BIG_INT_FOUR = Ref{BigInt}()

function __init__()
    BIG_INT_MINUS_ONE[] = BigInt(-1)
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
    # Add an inner constructor to overwrite the generic one that has underspecified type signature (for automatic conversion).
    # Automatic conversion would almost never work for Tv types, so it only causes ambiguities.
    PauliOperator{Tz,Tv}(phase, nqubits, xz::Tv) where {Tz<:AbstractArray{UInt8,0}, Tv<:AbstractVector{<:Unsigned}} = new(phase, nqubits, xz)
end

PauliOperator(phase::Tz, nqubits::Int, xz::Tv) where {Tz<:AbstractArray{UInt8,0}, Tv<:AbstractVector{<:Unsigned}} = PauliOperator{Tz,Tv}(phase, nqubits, xz)
PauliOperator(phase::UInt8, nqubits::Int, xz::Tv) where Tv<:AbstractVector{<:Unsigned} = PauliOperator(fill(UInt8(phase),()), nqubits, xz)
PauliOperator{Tz,Tv}(phase::UInt8, x::AbstractVector{Bool}, z::AbstractVector{Bool}) where {Tz, Tve<:Unsigned, Tv<:AbstractVector{Tve}} = PauliOperator(fill(UInt8(phase),()), length(x), vcat(reinterpret(Tve,BitVector(x).chunks),reinterpret(Tve,BitVector(z).chunks)))
PauliOperator(phase::UInt8, x::AbstractVector{Bool}, z::AbstractVector{Bool}) = PauliOperator{Array{UInt8,0},Vector{UInt}}(phase, x, z)
PauliOperator(x::AbstractVector{Bool}, z::AbstractVector{Bool}) = PauliOperator{Array{UInt8,0},Vector{UInt}}(0x0, x, z)
PauliOperator(xz::AbstractVector{Bool}) = PauliOperator{Array{UInt8,0},Vector{UInt}}(0x0, (@view xz[1:end÷2]), (@view xz[end÷2+1:end]))

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
    _P_str(a)
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

Base.show(io::IO, p::PauliOperator) = print(io, ["+ ","+i","- ","-i"][p.phase[]+1]*xz2str(xbit(p),zbit(p)))

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
# Stabilizers
##############################

abstract type AbstractQCState end # This could include classical bits

abstract type AbstractStabilizer <: AbstractQCState end # This includs only qubits

"""
Stabilizer, i.e. a list of commuting multi-qubit Hermitian Pauli operators.

Instances can be created with the `S` custom string macro or
as direct sum of other stabilizers.

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
struct Stabilizer{Tzv<:AbstractVector{UInt8}, Tm<:AbstractMatrix{<:Unsigned}} <: AbstractStabilizer
    phases::Tzv
    nqubits::Int
    xzs::Tm
end

function Stabilizer(paulis::AbstractVector{PauliOperator{Tz,Tv}}) where {Tz<:AbstractArray{UInt8,0},Tve<:Unsigned,Tv<:AbstractVector{Tve}}
    r = length(paulis)
    n = nqubits(paulis[1])
    stab = zero(Stabilizer{Vector{UInt8},Matrix{Tve}},r,n)
    for i in eachindex(paulis)
        stab[i] = paulis[i]
    end
    stab
end

Stabilizer(phases::AbstractVector{UInt8}, xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Stabilizer(
    phases, size(xs,2),
    vcat(hcat((BitArray(xs[i,:]).chunks for i in 1:size(xs,1))...)::Matrix{UInt},
         hcat((BitArray(zs[i,:]).chunks for i in 1:size(zs,1))...)::Matrix{UInt}) # type assertions to help Julia infer types
)

Stabilizer(phases::AbstractVector{UInt8}, xzs::AbstractMatrix{Bool}) = Stabilizer(phases, xzs[:,1:end÷2], xzs[:,end÷2+1:end])

Stabilizer(xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Stabilizer(zeros(UInt8, size(xs,1)), xs, zs)

Stabilizer(xzs::AbstractMatrix{Bool}) = Stabilizer(zeros(UInt8, size(xzs,1)), xzs[:,1:end÷2], xzs[:,end÷2+1:end])

Stabilizer(s::Stabilizer) = s

function _S_str(a)
    paulis = [_P_str(strip(s.match)) for s in eachmatch(r"[+-]?\h*[i]?\h*[XYZI_]+", a)]
    Stabilizer(paulis)
end

macro S_str(a)
    _S_str(a)
end

Base.getindex(stab::Stabilizer, i::Int) = PauliOperator(stab.phases[i], nqubits(stab), stab.xzs[:,i])
@inline function Base.getindex(stab::Stabilizer, r::Int, c::Int)
     # TODO this has code repetition with the Pauli getindex
     Tme = eltype(stab.xzs)
     x = (stab.xzs[_div(Tme,c-1)+1,r] & Tme(0x1)<<_mod(Tme,c-1))!=0x0
     z = (stab.xzs[end÷2+_div(Tme,c-1)+1,r] & Tme(0x1)<<_mod(Tme,c-1))!=0x0
     return (x,z)
end
Base.getindex(stab::Stabilizer, r) = Stabilizer(stab.phases[r], nqubits(stab), stab.xzs[:,r])
Base.getindex(stab::Stabilizer, r::Union{Colon,AbstractVector}, c::Union{Colon,AbstractVector}) = Stabilizer([s[c] for s in stab[r]])
Base.getindex(stab::Stabilizer, r::Union{Colon,AbstractVector}, c::Int) = stab[r,[c]]
Base.getindex(stab::Stabilizer, r::Int, c::Union{Colon,AbstractVector}) = stab[r][c]
Base.view(stab::Stabilizer, r) = Stabilizer(view(stab.phases, r), nqubits(stab), view(stab.xzs, :, r))

Base.iterate(stab::Stabilizer, state=1) = state>length(stab) ? nothing : (stab[state], state+1)

function Base.setindex!(stab::Stabilizer, pauli::PauliOperator, i)
    stab.phases[i] = pauli.phase[]
    #stab.xzs[:,i] = pauli.xz # TODO why is this assigment causing allocations
    for j in 1:length(pauli.xz)
        stab.xzs[j,i] = pauli.xz[j]
    end
    stab
end

function Base.setindex!(stab::Stabilizer, s::Stabilizer, i)
    stab.phases[i] = s.phases
    stab.xzs[:,i] = s.xzs
    stab
end

function Base.setindex!(stab::Stabilizer{Tzv,Tm}, (x,z)::Tuple{Bool,Bool}, i, j) where {Tzv<:AbstractVector{UInt8}, Tme<:Unsigned, Tm<:AbstractMatrix{Tme}} # TODO this has code repetitions with the Pauli setindex
    if x
        stab.xzs[_div(Tme,j-1)+1,        i] |= Tme(0x1)<<_mod(Tme,j-1)
    else
        stab.xzs[_div(Tme,j-1)+1,        i] &= ~(Tme(0x1)<<_mod(Tme,j-1))
    end
    if z
        stab.xzs[end÷2+_div(Tme,j-1)+1, i] |= Tme(0x1)<<_mod(Tme,j-1)
    else
        stab.xzs[end÷2+_div(Tme,j-1)+1, i] &= ~(Tme(0x1)<<_mod(Tme,j-1))
    end
    stab
end

Base.firstindex(stab::Stabilizer) = 1

Base.lastindex(stab::Stabilizer) = length(stab.phases)
Base.lastindex(stab::Stabilizer, i) = size(stab)[i]

Base.eachindex(stab::Stabilizer) = Base.OneTo(lastindex(stab.phases))

Base.axes(stab::Stabilizer) = (Base.OneTo(lastindex(stab)), Base.OneTo(nqubits(stab)))
Base.axes(stab::Stabilizer,i) = axes(stab)[i]

Base.size(stab::Stabilizer) = (length(stab.phases),nqubits(stab))
Base.size(stab::Stabilizer,i) = size(stab)[i]

Base.length(stab::Stabilizer) = length(stab.phases)

Base.show(io::IO, s::Stabilizer) = print(io,
                                         join([["+ ","+i","- ","-i"][s[i].phase[]+1]*xz2str(xbit(s[i]),zbit(s[i]))
                                               for i in eachindex(s)],
                                              '\n'))

Base.:(==)(l::Stabilizer, r::Stabilizer) = r.nqubits==l.nqubits && r.phases==l.phases && r.xzs==l.xzs

Base.hash(s::Stabilizer, h::UInt) = hash(s.nqubits, hash(s.phases, hash(s.xzs, h)))

Base.copy(s::Stabilizer) = Stabilizer(copy(s.phases), s.nqubits, copy(s.xzs))

function Base.zero(::Type{Stabilizer{Tzv, Tm}}, r, q) where {Tzv,T<:Unsigned,Tm<:AbstractMatrix{T}}
    phases = zeros(UInt8,r)
    xzs = zeros(UInt, _nchunks(q,T), r)
    Stabilizer(phases, q, xzs)::Stabilizer{Vector{UInt8},Matrix{UInt}}
end
Base.zero(::Type{S}, q) where {S<:Stabilizer}= zero(S, q, q)
Base.zero(::Type{Stabilizer}, r, q) = zero(Stabilizer{Vector{UInt8},Matrix{UInt}}, r, q)
Base.zero(s::S) where {S<:Stabilizer} = zero(S, length(s), nqubits(s))

"""Zero-out a given row of a [`Stabilizer`](@ref)"""
@inline function zero!(s::Stabilizer,i)
    s.xzs[:,i] .= zero(eltype(s.xzs))
    s.phases[i] = 0x0
    s
end

##############################
# Helpers for sublcasses of AbstractStabilizer that use Stabilizer as a tableau internally.
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

julia> tab(CliffordOperator(s))
+ Z
+ X

julia> typeof(tab(CliffordOperator(s)))
Stabilizer{Vector{UInt8}, Matrix{UInt64}}
```

See also: [`stabilizerview`](@ref), [`destabilizerview`](@ref), [`logicalxview`](@ref), [`logicalzview`](@ref)
"""
tab(s::Stabilizer) = s
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
struct Destabilizer{Tzv<:AbstractVector{UInt8},Tm<:AbstractMatrix{<:Unsigned}} <: AbstractStabilizer
    tab::Stabilizer{Tzv,Tm}
    function Destabilizer(s;noprocessing=false)
        r,n = size(s)
        if r==2n || noprocessing
            new{typeof(s.phases),typeof(s.xzs)}(s)
        elseif r<=n
            mixed_destab = MixedDestabilizer(s)
            tab = vcat(destabilizerview(mixed_destab),stabilizerview(mixed_destab))
            new{typeof(tab.phases),typeof(tab.xzs)}(tab)
        else
            throw(DomainError("To construct a `Destabilizer`, either use a 2n×n tableau of destabilizer and stabilizer rows that is directly used or an r×n (r<n) tableau of stabilizers from which destabilizers are automatically computed. Or better, just use `MixedDestabilizer`."))
        end
    end
end

function Base.show(io::IO, d::Destabilizer)
    show(io, destabilizerview(d))
    print(io, "\n━━" * "━"^size(d.tab,2) * "\n")
    show(io, stabilizerview(d))
end

Base.length(d::Destabilizer) = length(d.tab.phases)÷2

Base.copy(d::Destabilizer) = Destabilizer(copy(d.tab);noprocessing=true)

##############################
# Mixed Stabilizer states
##############################

"""
A slight improvement of the [`Stabilizer`](@ref) data structure that enables
more naturally and completely the treatment of mixed states, in particular when
the [`project!`](@ref) function is used.
"""
mutable struct MixedStabilizer{Tzv<:AbstractVector{UInt8}, Tm<:AbstractMatrix{<:Unsigned}} <: AbstractStabilizer
    tab::Stabilizer{Tzv,Tm} # TODO assert size on construction
    rank::Int
end

function MixedStabilizer(s::Stabilizer{Tzv,Tm}) where {Tzv<:AbstractVector{UInt8}, Tme<:Unsigned, Tm<:AbstractMatrix{Tme}}
    s, xr, zr = canonicalize!(s,ranks=true)
    spadded = zero(Stabilizer, nqubits(s))
    spadded[1:zr] = s[1:zr]
    MixedStabilizer(spadded,zr)
end

function Base.show(io::IO, ms::MixedStabilizer)
    println(io, "Rank $(ms.rank) stabilizer")
    show(io, stabilizerview(ms))
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
mutable struct MixedDestabilizer{Tzv<:AbstractVector{UInt8}, Tm<:AbstractMatrix{<:Unsigned}} <: AbstractStabilizer
    tab::Stabilizer{Tzv,Tm} # TODO assert size on construction
    rank::Int
end

# Added a lot of type assertions to help Julia infer types
function MixedDestabilizer(stab::Stabilizer{Tv,Tm}; undoperm=true) where {Tve,Tme,Tv<:AbstractVector{Tve},Tm<:AbstractMatrix{Tme}}
    rows,n = size(stab)
    stab, r, s, permx, permz = canonicalize_gott!(copy(stab))
    tab = zero(Stabilizer{Vector{Tve},Matrix{Tme}}, n*2, n)
    vstab = @view stab[1:r+s] # this view is necessary for cases of tableaux with redundant rows
    tab[n+1:n+r+s] = vstab # The Stabilizer part of the tableau
    for i in 1:r # The Destabilizer part
        tab[i,i] = (false,true)
    end
    for i in r+1:r+s
        tab[i,i] = (true,false)
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
        sX = Stabilizer(X)
        tab[r+s+1:n] = sX # One of the logical sets in the tableau
        A2 = H[1:r,r+s+1:end÷2]           #     r × n-r-s
        #Z = hcat(zeros(Bool,k,n),A2',zeros(Bool,k,s),I)
        Z = zeros(Bool, k, 2n)
        Z[:,n+1:n+r] .= A2'
        @inbounds @simd for i in 1:k Z[i,2n-k+i] = true end
        sZ = Stabilizer(Z)
        tab[n+r+s+1:end] = sZ # The other logical set in the tableau
    end
    if undoperm
        tab = tab[:,invperm(permx[permz])]::Stabilizer{Vector{Tve},Matrix{Tme}}
    end
    MixedDestabilizer(tab, r+s)::MixedDestabilizer{Vector{Tve},Matrix{Tme}}
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
    println(io, "Rank $(d.rank) stabilizer")
    show(io, destabilizerview(d))
    if d.rank != nqubits(d)
        print(io, "\n━━" * "━"^size(d.tab,2) * "\n")
        show(io, logicalxview(d))
        print(io, "\n━━" * "━"^size(d.tab,2) * "\n")
    else
        print(io, "\n══" * "═"^size(d.tab,2) * "\n")
    end
    show(io, stabilizerview(d))
    if d.rank != nqubits(d)
        print(io, "\n━━" * "━"^size(d.tab,2) * "\n")
        show(io, logicalzview(d))
    else
        print(io, "\n══" * "═"^size(d.tab,2) * "\n")
    end
end

Base.length(d::MixedDestabilizer) = length(d.tab)÷2

Base.copy(d::MixedDestabilizer) = MixedDestabilizer(copy(d.tab),d.rank)

##############################
# Subtableau views
##############################

"""A view of the subtableau corresponding to the stabilizer. See also [`tab`](@ref), [`destabilizerview`](@ref), [`logicalxview`](@ref), [`logicalzview`](@ref)"""
@inline stabilizerview(s::Stabilizer) = s
@inline stabilizerview(s::Destabilizer) = @view s.tab[end÷2+1:end]
@inline stabilizerview(s::MixedStabilizer) = @view s.tab[1:s.rank]
@inline stabilizerview(s::MixedDestabilizer) = @view s.tab[end÷2+1:end÷2+s.rank]

"""A view of the subtableau corresponding to the destabilizer. See also [`tab`](@ref), [`stabilizerview`](@ref), [`logicalxview`](@ref), [`logicalzview`](@ref)"""
@inline destabilizerview(s::Destabilizer) = @view s.tab[1:end÷2]
@inline destabilizerview(s::MixedDestabilizer) = @view s.tab[1:s.rank]

"""A view of the subtableau corresponding to the logical X operators. See also [`tab`](@ref), [`stabilizerview`](@ref), [`destabilizerview`](@ref), [`logicalzview`](@ref)"""
@inline logicalxview(s::MixedDestabilizer) = @view s.tab[s.rank+1:end÷2]
"""A view of the subtableau corresponding to the logical Z operators. See also [`tab`](@ref), [`stabilizerview`](@ref), [`destabilizerview`](@ref), [`logicalxview`](@ref)"""
@inline logicalzview(s::MixedDestabilizer) = @view s.tab[end÷2+s.rank+1:end]

"""The number of qubits of a given state."""
@inline nqubits(s::Stabilizer) = s.nqubits
@inline nqubits(s::AbstractStabilizer) = s.tab.nqubits

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

@inline function prodphase(l::PauliOperator, r::Stabilizer, i)::UInt8
    (l.phase[]+r.phases[i]+prodphase(l.xz, (@view r.xzs[:,i])))&0x3
end

@inline function prodphase(l::Stabilizer, r::PauliOperator, i)::UInt8
    (l.phases[i]+r.phase[]+prodphase((@view l.xzs[:,i]), r.xz))&0x3
end

@inline function prodphase(l::Stabilizer, r::Stabilizer, i, j)::UInt8
    (l.phases[i]+r.phases[j]+prodphase((@view l.xzs[:,i]), (@view r.xzs[:,j])))&0x3
end

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

@inline function comm(l::PauliOperator, r::Stabilizer, i::Int)::UInt8
    comm(l.xz,(@view r.xzs[:,i]))
end

@inline function comm(l::Stabilizer, r::PauliOperator, i::Int)::UInt8
    comm(r, l, i)
end

@inline function comm(s::Stabilizer, l::Int, r::Int)::UInt8
    comm((@view s.xzs[:,l]),(@view s.xzs[:,r]))
end

function comm(l::PauliOperator, r::Stabilizer)::Vector{UInt8}
    [comm(l,r,i) for i in 1:size(r,1)]
end

comm(l::Stabilizer, r::PauliOperator) = comm(r, l)

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

@inline function mul_left!(r::PauliOperator, l::Stabilizer, i::Int; phases::Val{B}=Val(true)) where B
    s = mul_left!(r.xz, (@view l.xzs[:,i]), phases=phases)
    B && (r.phase[] = (s+r.phase[]+l.phases[i])&0x3)
    r
end

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

@inline function rowswap!(s::Stabilizer, i, j; phases::Val{B}=Val(true)) where B # Written only so we can avoid copying in `canonicalize!`
    (i == j) && return
    B && begin s.phases[i], s.phases[j] = s.phases[j], s.phases[i] end
    @inbounds @simd for k in 1:size(s.xzs,1)
        s.xzs[k,i], s.xzs[k,j] = s.xzs[k,j], s.xzs[k,i]
    end
end

@inline function rowswap!(s::Destabilizer, i, j; phases::Val{B}=Val(true)) where B
    rowswap!(s.tab, i, j; phases)
    n = size(s.tab,1)÷2
    rowswap!(s.tab, i+n, j+n; phases)
end

@inline function rowswap!(s::MixedStabilizer, i, j; phases::Val{B}=Val(true)) where B
    rowswap!(s.tab, i, j; phases)
end

@inline function rowswap!(s::MixedDestabilizer, i, j; phases::Val{B}=Val(true)) where B
    rowswap!(s.tab, i, j; phases)
    n = nqubits(s)
    rowswap!(s.tab, i+n, j+n; phases)
end

@inline function mul_left!(s::Stabilizer, m, i; phases::Val{B}=Val(true)) where B
    extra_phase = mul_left!((@view s.xzs[:,m]), (@view s.xzs[:,i]); phases=phases)
    B && (s.phases[m] = (extra_phase+s.phases[m]+s.phases[i])&0x3)
    s
end

@inline function mul_left!(s::Destabilizer, i, j; phases::Val{B}=Val(true)) where B
    mul_left!(s.tab, j, i; phases=Val(false)) # Indices are flipped to preserve commutation constraints
    n = size(s.tab,1)÷2
    mul_left!(s.tab, i+n, j+n; phases=phases)
end

@inline function mul_left!(s::MixedStabilizer, i, j; phases::Val{B}=Val(true)) where B
    mul_left!(s.tab, i, j; phases=phases)
end

@inline function mul_left!(s::MixedDestabilizer, i, j; phases::Val{B}=Val(true)) where B
    mul_left!(s.tab, j, i; phases=Val(false)) # Indices are flipped to preserve commutation constraints
    n = nqubits(s)
    mul_left!(s.tab, i+n, j+n; phases=phases)
end

# TODO is this used anywhere?
"""Swap two columns of a stabilizer in place."""
@inline function colswap!(s::Stabilizer{Tzv,Tm}, i, j) where {Tzv<:AbstractVector{UInt8}, Tme<:Unsigned, Tm<:AbstractMatrix{Tme}}
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

function colpermute!(s::Stabilizer, perm) # TODO rename and make public, same as permute and maybe Base.permute!
    for r in 1:size(s,1)
        s[r] = s[r][perm] # TODO make a local temporary buffer row instead of constantly allocating new rows
    end
    s
end

function ishermitian() # TODO write it both for paulis and stabilizers (ugh... stabilizers are always so)
end

function check_allrowscommute(stabilizer::Stabilizer)
    for i in eachindex(stabilizer)
        for j in eachindex(stabilizer)
            i==j && continue
            comm(stabilizer[i],stabilizer[j])==0x0 || return false
        end
    end
    return true
end

function Base.vcat(stabs::Stabilizer...)
    Stabilizer(vcat((s.phases for s in stabs)...),
               stabs[1].nqubits,
               hcat((s.xzs for s in stabs)...))
end

##############################
# Unitary Clifford Operations
##############################

function Base.:(*)(p::AbstractCliffordOperator, s::AbstractStabilizer; phases::Bool=true)
    s = copy(s)
    _apply!(s,p; phases=Val(phases))
end
function apply!(stab::AbstractStabilizer, op::AbstractCliffordOperator; phases::Bool=true)
    _apply!(stab,op; phases=Val(phases))
end
function apply!(stab::AbstractStabilizer, op::AbstractCliffordOperator, indices; phases::Bool=true)
    _apply!(stab,op,indices; phases=Val(phases))
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

"""The F(2,2) matrix of a given stabilizer, represented as the concatenation of two binary matrices, one for X and one for Z."""
function stab_to_gf2(s::Stabilizer)::Matrix{Bool}
    r, n = size(s)
    H = zeros(Bool,r,2n)
    for iᵣ in 1:r
        @inbounds @simd for iₙ in 1:n
            H[iᵣ,iₙ], H[iᵣ,iₙ+n] = s[iᵣ,iₙ]
        end
    end
    H
end

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
    all((ut[i, i] for i in 1:size(ut,1)))
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
        i ∈ goodindices && continue
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
    c = canonicalize!(copy(s))
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
    c = canonicalize!(copy(s))
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
include("misc_ops.jl")
include("classical_register.jl")
include("enumeration.jl")
include("randoms.jl")
include("useful_states.jl")
include("experimental/Experimental.jl")
include("graphs.jl")
include("entanglement.jl")
include("precompiles.jl")

end #module
