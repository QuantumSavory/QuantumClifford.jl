"""
A module for simulation of Clifford circuits.
"""
module SimpleClifford

# TODO document phases=false

# TODO PauliOperator should be mutable so that the phase is not a zero-dim
# array, but to do that without sacrificing functionality we need to implement
# views. That would also fix the current issue of getindex returning a view,
# instead of us having a dedicated view type.

# TODO remove most of the convenience properties and make helper functions.

# TODO Operations between Clifford operators are very slow

import LinearAlgebra
import Random
import RecipesBase

export @P_str, PauliOperator, ⊗, I, X, Y, Z, permute,
    @S_str, Stabilizer, prodphase, comm, ⊕, check_allrowscommute,
    Destabilizer, MixedStabilizer, MixedDestabilizer,
    nqubits, stabilizerview, destabilizerview, logicalxview, logicalzview,
    canonicalize!, canonicalize_rref!, canonicalize_gott!, colpermute!,
    generate!, project!, reset_qubits!, traceout!,
    apply!,
    CliffordOperator, @C_str, CNOT, SWAP, Hadamard, Phase, CliffordId,
    tensor_pow,
    stab_to_gf2, gf2_gausselim!, gf2_isinvertible, gf2_invert, gf2_H_to_G,
    perm_inverse, perm_product,
    single_z, single_x,
    random_invertible_gf2,
    random_pauli, random_stabilizer, random_singlequbitop,
    BadDataStructure

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

abstract type AbstractCliffordOperator end

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
in a single concatenated padded array of UInt64 chunks of a bit array.

```jldoctest
julia> p = P"-IZXY";

julia> p.xz
2-element Array{UInt64,1}:
 0x000000000000000c
 0x000000000000000a
```
"""
struct PauliOperator{Tz<:AbstractArray{UInt8,0},Tv<:AbstractVector{UInt64}} <: AbstractCliffordOperator
    phase::Tz
    nqubits::Int
    xz::Tv
end

PauliOperator(phase::UInt8, nqubits::Int, xz::Tv) where Tv<:AbstractVector{UInt64} = PauliOperator(fill(phase,()), nqubits, xz)
PauliOperator(phase::UInt8, x::T, z::T) where T<:AbstractVector{Bool} = PauliOperator(fill(phase,()), length(x), vcat(BitVector(x).chunks,BitVector(z).chunks))

function Base.getproperty(p::PauliOperator, name::Symbol)
    if name==:xview
        @view p.xz[1:end÷2]
    elseif name==:zview
        @view p.xz[end÷2+1:end]
    elseif name==:xbit
        b = BitArray(UndefInitializer(),(p.nqubits,))
        b.chunks = p.xview
        b
    elseif name==:zbit
        b = BitArray(UndefInitializer(),(p.nqubits,))
        b.chunks = p.zview
        b
    else
        getfield(p, name)
    end
end

Base.propertynames(p::PauliOperator, private=false) = (:phase,:nqubits,:xz,:xbit,:zbit,:xview,:zview)

macro P_str(a)
    letters = filter(x->occursin(x,"_IZXY"),a)
    phase = phasedict[strip(filter(x->!occursin(x,"_IZXY"),a))]
    PauliOperator(phase, [l=='X'||l=='Y' for l in letters], [l=='Z'||l=='Y' for l in letters])
end

Base.getindex(p::PauliOperator, i::Int) = (p.xz[_div64(i-1)+1] & UInt64(0x1)<<_mod64(i-1))!=0x0, (p.xz[end>>1+_div64(i-1)+1] & UInt64(0x1)<<_mod64(i-1))!=0x0
Base.getindex(p::PauliOperator, r) = PauliOperator(p.phase[], p.xbit[r], p.zbit[r])

function Base.setindex!(p::PauliOperator, (x,z)::Tuple{Bool,Bool}, i)
    if x
        p.xz[_div64(i-1)+1] |= UInt64(0x1)<<_mod64(i-1)
    else
        p.xz[_div64(i-1)+1] &= ~(UInt64(0x1)<<_mod64(i-1))
    end
    if z
        p.xz[end>>1+_div64(i-1)+1] |= UInt64(0x1)<<_mod64(i-1)
    else
        p.xz[end>>1+_div64(i-1)+1] &= ~(UInt64(0x1)<<_mod64(i-1))
    end
    p
end

Base.firstindex(p::PauliOperator) = 1

Base.lastindex(p::PauliOperator) = p.nqubits

Base.eachindex(p::PauliOperator) = 1:p.nqubits

Base.size(pauli::PauliOperator) = (pauli.nqubits,)

Base.length(pauli::PauliOperator) = pauli.nqubits

xz2str(x,z) = join(toletter[e] for e in zip(x,z))

Base.show(io::IO, p::PauliOperator) = print(io, ["+ ","+i","- ","-i"][p.phase[]+1]*xz2str(p.xbit,p.zbit))

Base.:(==)(l::PauliOperator, r::PauliOperator) = r.phase==l.phase && r.nqubits==l.nqubits && r.xz==l.xz

Base.hash(p::PauliOperator, h::UInt) = hash((p.phase,p.nqubits,p.xz), h)

Base.copy(p::PauliOperator) = PauliOperator(copy(p.phase),p.nqubits,copy(p.xz))

##############################
# Stabilizers
##############################

abstract type AbstractStabilizer end

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

julia> s⊕s
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

There are no automatic checks for correctness (i.e. independence of all rows,
commutativity of all rows, hermiticity of all rows). The rank (number of rows)
is permitted to be less than the number of qubits (number of columns):
canonilization, projection, etc. continue working in that case. To great extent
this library uses the `Stabilizer` data structure simply as a tableau. This
might be properly abstracted away in future versions.

See also: [`PauliOperator`](@ref), [`canonicalize!`](@ref)
"""
struct Stabilizer{Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{UInt64}} <: AbstractStabilizer
    phases::Tv
    nqubits::Int
    xzs::Tm
end

Stabilizer(paulis::AbstractVector{PauliOperator{Tz,Tv}}) where {Tz<:AbstractArray{UInt8,0},Tv<:AbstractVector{UInt64}} = Stabilizer(vcat((p.phase for p in paulis)...), paulis[1].nqubits, vcat((p.xz' for p in paulis)...))

Stabilizer(phases::AbstractVector{UInt8}, xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Stabilizer(
    phases, size(xs,2),
    hcat(vcat((BitArray(xs[i,:]).chunks' for i in 1:size(xs,1))...),
         vcat((BitArray(zs[i,:]).chunks' for i in 1:size(zs,1))...))
)

Stabilizer(phases::AbstractVector{UInt8}, xzs::AbstractMatrix{Bool}) = Stabilizer(phases, xzs[:,1:end÷2], xzs[:,end÷2+1:end])

Stabilizer(xs::AbstractMatrix{Bool}, zs::AbstractMatrix{Bool}) = Stabilizer(zeros(UInt8, size(xs,1)), xs, zs)

Stabilizer(xzs::AbstractMatrix{Bool}) = Stabilizer(zeros(UInt8, size(xzs,1)), xzs[:,1:end÷2], xzs[:,end÷2+1:end])

macro S_str(a)
    paulis = [eval(quote @P_str($(strip(s))) end) for s in split(a,'\n')] #TODO seriously!?
    Stabilizer(paulis)
end

Base.getindex(stab::Stabilizer, i::Int) = PauliOperator((@view stab.phases[i]), nqubits(stab), (@view stab.xzs[i,:]))
Base.getindex(stab::Stabilizer, r) = Stabilizer((@view stab.phases[r]), nqubits(stab), (@view stab.xzs[r,:]))
Base.getindex(stab::Stabilizer, r, c) = Stabilizer([s[c] for s in stab[r]])

Base.iterate(stab::Stabilizer, state=1) = state>length(stab) ? nothing : (stab[state], state+1)

function Base.setindex!(stab::Stabilizer, pauli::PauliOperator, i)
    stab.phases[i] = pauli.phase[]
    stab.xzs[i,:] = pauli.xz
    stab
end

function Base.setindex!(stab::Stabilizer, s::Stabilizer, i)
    stab.phases[i] = s.phases
    stab.xzs[i,:] = s.xzs
    stab
end

function Base.setindex!(stab::Stabilizer, (x,z)::Tuple{Bool,Bool}, i, j)
    if x
        stab.xzs[i,_div64(j-1)+1] |= UInt64(0x1)<<_mod64(j-1)
    else
        stab.xzs[i,_div64(j-1)+1] &= ~(UInt64(0x1)<<_mod64(j-1))
    end
    if z
        stab.xzs[i,end>>1+_div64(j-1)+1] |= UInt64(0x1)<<_mod64(j-1)
    else
        stab.xzs[i,end>>1+_div64(j-1)+1] &= ~(UInt64(0x1)<<_mod64(j-1))
    end
    stab
end

Base.firstindex(stab::Stabilizer) = 1

Base.lastindex(stab::Stabilizer) = length(stab.phases)

Base.eachindex(stab::Stabilizer) = 1:length(stab.phases)

Base.size(stab::Stabilizer) = (length(stab.phases),nqubits(stab))
Base.size(stab::Stabilizer,i) = size(stab)[i]

Base.length(stab::Stabilizer) = length(stab.phases)

Base.show(io::IO, s::Stabilizer) = print(io,
                                         join([["+ ","+i","- ","-i"][s[i].phase[]+1]*xz2str(s[i].xbit,s[i].zbit)
                                               for i in eachindex(s)],
                                              '\n'))

Base.:(==)(l::Stabilizer, r::Stabilizer) = r.nqubits==l.nqubits && r.phases==l.phases && r.xzs==l.xzs

Base.hash(s::Stabilizer, h::UInt) = hash(s.nqubits, s.phases, s.xzs, h)

Base.copy(s::Stabilizer) = Stabilizer(copy(s.phases), s.nqubits, copy(s.xzs))

##############################
# Helpers for sublcasses of AbstractStabilizer that use Stabilizer as a tableau internally.
##############################

Base.:(==)(l::T, r::S) where {T<:AbstractStabilizer, S<:AbstractStabilizer} = T==S && r.tab==l.tab

Base.hash(s::T, h::UInt) where {T<:AbstractStabilizer} = hash(T, s.tab, h)

function apply!(s::AbstractStabilizer, p::AbstractCliffordOperator; phases::Bool=true)
    apply!(s.tab,p; phases=phases)
    s
end

##############################
# Destabilizer formalism
##############################

function destabilizer_generators(stab::Stabilizer)::Tuple{Stabilizer,Stabilizer}
    stab = canonicalize!(copy(stab))
    dest = zero(stab)
    s = 1
    e = size(stab.xzs,2)>>1
    op = (false, true)
    for i in eachindex(stab)
        j = unsafe_bitfindnext_(stab.xzs[i,s:e],1)
        if isnothing(j)
            s = e+1
            e = 2*e
            op = (true, false)
            j = unsafe_bitfindnext_(stab.xzs[i,s:e],1)
        end
        dest[i,j] = op
    end
    dest, stab
end

"""
A tableau representation of a pure stabilizer state. The tableau tracks the
destabilizers as well, for efficient projections. On initialization there are
no checks that the provided state is indeed pure. This enables the use of this
data structure for mixed stabilizer state, but a better choice would be to use
[`MixedDestabilizer`](@ref).
"""
struct Destabilizer{Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{UInt64}} <: AbstractStabilizer
    tab::Stabilizer{Tv,Tm}
    function Destabilizer(s;noprocessing=false)
        if noprocessing
            new{typeof(s.phases),typeof(s.xzs)}(s)
        else
            tab = vcat(destabilizer_generators(s)...)
            new{typeof(s.phases),typeof(s.xzs)}(tab)
        end
    end
end

function Base.show(io::IO, d::Destabilizer)
    show(io, destabilizerview(d))
    print(io, "\n━━" * "━"^size(d.tab,2) * "\n")
    show(io, stabilizerview(d))
end

Base.copy(d::Destabilizer) = Destabilizer(copy(d.tab);noprocessing=true)

##############################
# Mixed Stabilizer states
##############################

"""
A slight improvement of the [`Stabilizer`](@ref) data structure that enables
more naturally and completely the treatment of mixed states, in particular when
the [`project!`](@ref) function is used.
"""
mutable struct MixedStabilizer{Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{UInt64}} <: AbstractStabilizer
    tab::Stabilizer{Tv,Tm} # TODO assert size on construction
    rank::Int
end

function MixedStabilizer(s::Stabilizer)
    s = canonicalize!(s)
    rp1 = findfirst(mapslices(row->all(==(UInt64(0)),row),s.xzs; dims=(2,)))
    r = isnothing(rp1) ? size(s, 1) : rp1-1
    spadded = zero(Stabilizer, nqubits(s))
    spadded[1:r] = s
    MixedStabilizer(spadded,r)
end

function Base.show(io::IO, ms::MixedStabilizer)
    println(io, "Rank $(ms.rank) stabilizer")
    show(io, stabilizerview(ms))
end

Base.copy(ms::MixedStabilizer) = MixedStabilizer(copy(ms.tab),ms.rank)

##############################
# Mixed Destabilizer states
##############################

"""
A tableau representation for mixed stabilizer states that keeps track of the
destabilizers in order to provide efficient projection operations.
"""
mutable struct MixedDestabilizer{Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{UInt64}} <: AbstractStabilizer
    tab::Stabilizer{Tv,Tm} # TODO assert size on construction
    rank::Int
end

function MixedDestabilizer(stab::Stabilizer; undoperm=true)
    r,n = size(stab)
    r==n && return MixedDestabilizer(Destabilizer(stab).tab, r)
    stab, r, s, permx, permz = canonicalize_gott!(stab)
    n = nqubits(stab)
    tab = zero(Stabilizer, n*2, n)
    tab[n+1:n+r+s] = stab
    for i in 1:r
        tab[i,i] = (false,true)
    end
    for i in r+1:r+s
        tab[i,i] = (true,false)
    end
    H = stab_to_gf2(stab)
    k = n - r - s
    E = H[r+1:end,end÷2+r+s+1:end]
    C1 = H[1:r,end÷2+r+1:end÷2+r+s]
    C2 = H[1:r,end÷2+r+s+1:end]
    i = LinearAlgebra.I
    U2 = E'
    V1 = (E' * C1' + C2').%2 .!= 0x0
    X = hcat(zeros(Bool,k,r),U2,i,V1,zeros(Bool,k,s+k))
    sX = Stabilizer(X)
    tab[r+s+1:n] = sX
    A2 = H[1:r,r+s+1:end÷2]
    Z = hcat(zeros(Bool,k,n),A2',zeros(Bool,k,s),i)
    sZ = Stabilizer(Z)
    tab[n+r+s+1:end] = sZ
    if undoperm
        tab = tab[:,perm_inverse(permx[permz])]
    end
    MixedDestabilizer(tab, r+s)
end

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

Base.copy(d::MixedDestabilizer) = MixedDestabilizer(copy(d.tab),d.rank)

##############################
# Subtableau views
##############################

@inline stabilizerview(s::Stabilizer) = s
@inline stabilizerview(s::Destabilizer) = s.tab[end÷2+1:end]
@inline stabilizerview(s::MixedStabilizer) = s.tab[1:s.rank]
@inline stabilizerview(s::MixedDestabilizer) = s.tab[end÷2+1:end÷2+s.rank]

@inline destabilizerview(s::Destabilizer) = s.tab[1:end÷2]
@inline destabilizerview(s::MixedDestabilizer) = s.tab[1:s.rank]

@inline logicalxview(s::MixedDestabilizer) = s.tab[s.rank+1:end÷2]
@inline logicalzview(s::MixedDestabilizer) = s.tab[end÷2+s.rank+1:end]

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
@inline function prodphase(l::AbstractVector{UInt64}, r::AbstractVector{UInt64})::UInt64
    res = 0
    len = length(l)>>1
    @inbounds @simd for i in 1:len
        x1, x2, z1, z2 = l[i], r[i], l[i+len], r[i+len]
        res += count_ones((~z2 & x2 & ~x1 & z1) | ( z2 & ~x2 & x1 & z1) | (z2 &  x2 & x1 & ~z1))
        res -= count_ones(( z2 & x2 & ~x1 & z1) | (~z2 &  x2 & x1 & z1) | (z2 & ~x2 & x1 & ~z1))
    end
    unsigned(res)
end

@inline function prodphase(l::PauliOperator, r::PauliOperator)::UInt8
    (l.phase[]+r.phase[]+prodphase(l.xz,r.xz))&0x3
end

@inline function prodphase(l::PauliOperator, r::Stabilizer, i)::UInt8 # TODO rewrite it in a way that does not use views in order to have fewer allocations
    (l.phase[]+r.phases[i]+prodphase(l.xz, (@view r.xzs[i,:])))&0x3
end

@inline function prodphase(l::Stabilizer, r::PauliOperator, i)::UInt8
    (l.phases[i]+r.phase[]+prodphase((@view l.xzs[i,:]), r.xz))&0x3
end

@inline function prodphase(l::Stabilizer, r::Stabilizer, i, j)::UInt8
    (l.phases[i]+r.phases[j]+prodphase((@view l.xzs[i,:]), (@view r.xzs[j,:])))&0x3
end

@inline function xor_bits_(v::UInt64)
    v ⊻= v >> 32
    v ⊻= v >> 16
    v ⊻= v >> 8
    v ⊻= v >> 4
    v ⊻= v >> 2
    v ⊻= v >> 1
    return v&1
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
@inline function comm(l::AbstractVector{UInt64}, r::AbstractVector{UInt64})::UInt8
    res = UInt64(0)
    len = length(l)>>1
    @inbounds @simd for i in 1:len
        res ⊻= (l[i+len] & r[i]) ⊻ (l[i] & r[i+len])
    end
    xor_bits_(res)
end

@inline function comm(l::PauliOperator, r::PauliOperator)::UInt8
    comm(l.xz,r.xz)
end

@inline function comm(l::PauliOperator, r::Stabilizer, i::Int)::UInt8
    comm(l.xz,(@view r.xzs[i,:]))
end

function comm(l::PauliOperator, r::Stabilizer)::Vector{UInt8}
    [comm(l,r,i) for i in 1:size(r,1)]
end

function Base.:(*)(l::PauliOperator, r::PauliOperator)
    PauliOperator(prodphase(l,r), l.nqubits, l.xz .⊻ r.xz)
end

(⊗)(l::PauliOperator, r::PauliOperator) = PauliOperator((l.phase[]+r.phase[])&0x3, vcat(l.xbit,r.xbit), vcat(l.zbit,r.zbit))

function Base.:(*)(l, r::PauliOperator)
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

@inline function rowswap!(s::Stabilizer, i, j; phases::Bool=true) # Written only so we can avoid copying in `canonicalize!`
    (i == j) && return
    phases && begin s.phases[i], s.phases[j] = s.phases[j], s.phases[i] end
    @inbounds @simd for k in 1:size(s.xzs,2)
        s.xzs[i,k], s.xzs[j,k] = s.xzs[j,k], s.xzs[i,k]
    end
end

@inline function rowswap!(s::Destabilizer, i, j; phases::Bool=true)
    rowswap!(s.tab, i, j; phases=phases)
    n = size(s.tab,1)÷2
    rowswap!(s.tab, i+n, j+n; phases=phases)
end

@inline function rowswap!(s::MixedStabilizer, i, j; phases::Bool=true)
    rowswap!(s.tab, i, j; phases=phases)
end

@inline function rowswap!(s::MixedDestabilizer, i, j; phases::Bool=true)
    rowswap!(s.tab, i, j; phases=phases)
    n = nqubits(s)
    rowswap!(s.tab, i+n, j+n; phases=phases)
end

@inline function rowmul!(s::Stabilizer, m, i; phases::Bool=true)
    @inbounds @simd for d in 1:size(s.xzs,2)
        s.xzs[m,d] ⊻= s.xzs[i,d]
    end
    phases && (s.phases[m] = prodphase(s,s,m,i))
end

@inline function rowmul!(s::Destabilizer, i, j; phases::Bool=true)
    rowmul!(s.tab, j, i; phases=false)
    n = size(s.tab,1)÷2
    rowmul!(s.tab, i+n, j+n; phases=phases)
end

@inline function rowmul!(s::MixedStabilizer, i, j; phases::Bool=true)
    rowmul!(s.tab, i, j; phases=phases)
end

@inline function rowmul!(s::MixedDestabilizer, i, j; phases::Bool=true)
    rowmul!(s.tab, j, i; phases=false)
    n = nqubits(s)
    rowmul!(s.tab, i+n, j+n; phases=phases)
end

# TODO document as a more efficient way of swaping only two columns (instead of using permutation indexing)
@inline function colswap!(s::Stabilizer, i, j)
    lowbit = UInt64(1)
    ibig = _div64(i-1)+1
    ismall = _mod64(i-1)
    ismallm = lowbit<<(ismall)
    jbig = _div64(j-1)+1
    jsmall = _mod64(j-1)
    jsmallm = lowbit<<(jsmall)
    for off in [0,size(s.xzs,2)÷2]
        ibig += off
        jbig += off
        @inbounds for k in 1:size(s.xzs,1)
            ival = s.xzs[k,ibig] & ismallm
            jval = s.xzs[k,jbig] & jsmallm
            s.xzs[k,ibig] &= ~ismallm
            s.xzs[k,jbig] &= ~jsmallm
            if ismall>jsmall
                s.xzs[k,ibig] |= jval<<(ismall-jsmall)
                s.xzs[k,jbig] |= ival>>(ismall-jsmall)
            elseif ismall<jsmall
                s.xzs[k,ibig] |= jval>>(jsmall-ismall)
                s.xzs[k,jbig] |= ival<<(jsmall-ismall)
            else
                s.xzs[k,ibig] |= jval
                s.xzs[k,jbig] |= ival
            end
        end
    end
end

# Copied from base/bitarray.jl
const _msk64 = ~UInt64(0)
@inline _div64(l) = l >> 6
@inline _mod64(l) = l & 63
function unsafe_bitfindnext_(chunks::AbstractVector{UInt64}, start::Integer)
    chunk_start = _div64(start-1)+1
    within_chunk_start = _mod64(start-1)
    mask = _msk64 << within_chunk_start

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

"""
Canonicalize a stabilizer (in place).

Assumes the input is a valid stabilizer (all operators commute and have
real phases). It permits redundant generators and identity generators.

```jldoctest
julia> ghz = S"XXXX
               ZZII
               IZZI
               IIZZ";

julia> canonicalize!(ghz)
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ

julia> canonicalize!(S"XXXX
                       IZZI
                       IIZZ")
+ XXXX
+ _Z_Z
+ __ZZ

julia> canonicalize!(S"XXXX
                       ZZII
                       IZZI
                       IZIZ
                       IIZZ")
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ
+ ____
```

Based on arxiv:1210.6646.
See arxiv:0505036 for other types of canonicalization.
"""
function canonicalize!(stabilizer::Stabilizer; phases::Bool=true)
    xzs = stabilizer.xzs
    xs = @view xzs[:,1:end÷2]
    zs = @view xzs[:,end÷2+1:end]
    lowbit = UInt64(0x1)
    zero64 = UInt64(0x0)
    rows, columns = size(stabilizer)
    i = 1
    for j in 1:columns
        # find first row with X or Y in col `j`
        jbig = _div64(j-1)+1
        jsmall = lowbit<<_mod64(j-1)
        k = findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                      (@view xs[i:end,jbig]))
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases=phases)
            for m in 1:rows
                if xs[m,jbig]&jsmall!=zero64 && m!=i # if X or Y
                    rowmul!(stabilizer, m, i; phases=phases)
                end
            end
            i += 1
        end
    end
    for j in 1:columns
        # find first row with Z in col `j`
        jbig = _div64(j-1)+1
        jsmall = lowbit<<_mod64(j-1)
        k = findfirst(e->e&(jsmall)!=zero64,
                      (@view zs[i:end,jbig]))
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases=phases)
            for m in 1:rows
                if zs[m,jbig]&jsmall!=zero64 && m!=i # if Z or Y
                    rowmul!(stabilizer, m, i; phases=phases)
                end
            end
            i += 1
        end
    end
    stabilizer
end

function canonicalize!(ms::MixedStabilizer; phases::Bool=true)
    canonicalize!(stabilizerview(ms); phases=phases)
end

"""
Canonicalize a stabilizer (in place) along only some columns.

This uses different canonical form from `canonicalize!`. It also indexes in
reverse in order to make its use in `traceout!` more efficient.

Based on arxiv:0505036.
"""
function canonicalize_rref!(state::AbstractStabilizer, colindices::AbstractVector{T}; phases::Bool=true) where {T<:Integer}
    xzs = stabilizerview(state).xzs
    xs = @view xzs[:,1:end÷2]
    zs = @view xzs[:,end÷2+1:end]
    lowbit = UInt64(0x1)
    zero64 = UInt64(0x0)
    rows, columns = size(stabilizerview(state))
    i = rows
    for j in colindices
        jbig = _div64(j-1)+1
        jsmall = lowbit<<_mod64(j-1)
        k = findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                      (@view xs[1:i,jbig]))
        if k !== nothing
            rowswap!(state, k, i; phases=phases)
            for m in 1:rows
                if xs[m,jbig]&jsmall!=zero64 && m!=i # if X or Y
                    rowmul!(state, m, i; phases=phases)
                end
            end
            i -= 1
        end
        k = findfirst(e->e&(jsmall)!=zero64,
                      (@view zs[1:i,jbig]))
        if k !== nothing
            rowswap!(state, k, i; phases=phases)
            for m in 1:rows
                if zs[m,jbig]&jsmall!=zero64 && m!=i # if Z or Y
                    rowmul!(state, m, i; phases=phases)
                end
            end
            i -= 1
        end
    end
    state, i
end

function gott_standard_form_indices(chunks2D, rows, cols; skip=0)
    goodindices = Int[]
    j = 1
    r = 1
    for r in skip+1:rows
        i = unsafe_bitfindnext_(chunks2D[r,:],skip+1)
        isnothing(i) && break
        i ∈ goodindices && continue
        push!(goodindices, i)
    end
    rank = length(goodindices)
    if rank>0
        badindices = [r for r in 1+skip:goodindices[end] if !(r ∈ goodindices)]
        return vcat(1:skip, goodindices, badindices, goodindices[end]+1:cols), rank
    else
        return 1:cols, rank
    end
end

function colpermute!(s::Stabilizer, perm) # TODO rename and make public, same as permute and maybe Base.permute!
    for r in 1:size(s,1)
        s[r] = s[r][perm]
    end
    s
end

function canonicalize_gott!(stabilizer::Stabilizer; phases::Bool=true)
    xzs = stabilizer.xzs
    xs = @view xzs[:,1:end÷2]
    zs = @view xzs[:,end÷2+1:end]
    lowbit = UInt64(0x1)
    zero64 = UInt64(0x0)
    rows, columns = size(stabilizer)
    i = 1
    for j in 1:columns
        # find first row with X or Y in col `j`
        jbig = _div64(j-1)+1
        jsmall = lowbit<<_mod64(j-1)
        k = findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                      xs[i:end,jbig])
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases=phases)
            for m in 1:rows
                if xs[m,jbig]&jsmall!=zero64 && m!=i # if X or Y
                    rowmul!(stabilizer, m, i; phases=phases)
                end
            end
            i += 1
        end
    end
    xperm, r = gott_standard_form_indices((@view xzs[:,1:end÷2]),rows,columns)
    colpermute!(stabilizer,xperm)
    i = r+1
    for j in r+1:columns
        # find first row with Z in col `j`
        jbig = _div64(j-1)+1
        jsmall = lowbit<<_mod64(j-1)
        k = findfirst(e->e&(jsmall)!=zero64,
                      zs[i:end,jbig])
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases=phases)
            for m in 1:rows
                if zs[m,jbig]&jsmall!=zero64 && m!=i # if Z or Y
                    rowmul!(stabilizer, m, i; phases=phases)
                end
            end
            i += 1
        end
    end
    zperm, s = gott_standard_form_indices((@view xzs[:,end÷2+1:end]),rows,columns,skip=r)
    colpermute!(stabilizer,zperm)
    stabilizer, r, s, xperm, zperm
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

function ⊕(l::Stabilizer, r::Stabilizer)
    lone = zero(l[1])
    rone = zero(r[1])
    paulis = vcat([l[i]⊗rone for i in eachindex(l)],
                  [lone⊗r[i] for i in eachindex(r)]
                 )
    Stabilizer(paulis)
end

function Base.vcat(stabs::Stabilizer...)
    Stabilizer(vcat((s.phases for s in stabs)...),
               stabs[1].nqubits,
               vcat((s.xzs for s in stabs)...))
end

##############################
# Projections and Measurements
##############################

"""
Generate a Pauli operator by using operators from a given the Stabilizer.

**It assumes the stabilizer is already canonicalized.** It modifies
the Pauli operator in place. It assumes the operator can be generated up to a phase.
That phase is left in the modified operator, which should be the identity up to a phase.
Returns the new operator and the list of indices denoting the elements of
`stabilizer` that were used for the generation.

```jldoctest
julia> ghz = S"XXXX
               ZZII
               IZZI
               IIZZ";

julia> canonicalize!(ghz)
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ

julia> generate!(P"-ZIZI", ghz)
(- ____, [2, 4])
```
"""
function generate!(pauli::PauliOperator, stabilizer::Stabilizer; phases::Bool=true, saveindices::Bool=true) # TODO there is stuff that can be abstracted away here and in canonicalize!
    rows, columns = size(stabilizer)
    xzs = stabilizer.xzs
    xs = @view xzs[:,1:end÷2]
    zs = @view xzs[:,end÷2+1:end]
    lowbit = UInt64(0x1)
    zero64 = UInt64(0x0)
    px,pz = pauli.xview, pauli.zview
    used_indices = Int[]
    used = 0
    # remove Xs
    while (i=unsafe_bitfindnext_(px,1)) !== nothing
        jbig = _div64(i-1)+1
        jsmall = lowbit<<_mod64(i-1)
        candidate = findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                              xs[used+1:end,jbig])
        if isnothing(candidate)
            return nothing
        else
            used += candidate
        end
        # TODO, this is just a long explicit way to write it... learn more about broadcast
        phases && (pauli.phase[] = prodphase(pauli,stabilizer,used))
        pauli.xz .⊻= @view xzs[used,:]
        saveindices && push!(used_indices, used)
    end
    # remove Zs
    while (i=unsafe_bitfindnext_(pz,1)) !== nothing
        jbig = _div64(i-1)+1
        jsmall = lowbit<<_mod64(i-1)
        candidate = findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                              zs[used+1:end,jbig])
        if isnothing(candidate)
            return nothing
        else
            used += candidate
        end
        # TODO, this is just a long explicit way to write it... learn more about broadcast
        phases && (pauli.phase[] = prodphase(pauli,stabilizer,used))
        pauli.xz .⊻= @view xzs[used,:]
        saveindices && push!(used_indices, used)
    end
    pauli, used_indices
end

"""
Project the state of a Stabilizer on the two eigenspaces of a Pauli operator.

Assumes the input is a valid stabilizer.
The projection is done inplace on that stabilizer and it does not modify the
projection operator.

It returns

 - a stabilizer that might not be in canonical form
 - the index of the row where the non-commuting operator was (that row is now equal to `pauli`; its phase is not updated and for a faithful measurement simulation it needs to be randomized by the user)
 - and the result of the projection if there was no non-cummuting operator (`nothing` otherwise)

If `keep_result==false` that result of the projection in case of anticommutation
is not computed, sparing a canonicalization operation.

Here is an example of a projection destroing entanglement:

```jldoctest
julia> ghz = S"XXXX
               ZZII
               IZZI
               IIZZ";

julia> canonicalize!(ghz)
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ

julia> state, anticom_index, result = project!(ghz, P"ZIII");

julia> state
+ Z___
+ Z__Z
+ _Z_Z
+ __ZZ

julia> canonicalize!(state)
+ Z___
+ _Z__
+ __Z_
+ ___Z

julia> anticom_index, result
(1, nothing)
```

And an example of projection consistent with the stabilizer state.

```jldoctest
julia> s = S"ZII
             IXI
             IIY";

julia> canonicalize!(s)
+ _X_
+ __Y
+ Z__

julia> state, anticom_index, result = project!(s, P"-ZII");

julia> state
+ _X_
+ __Y
+ Z__

julia> anticom_index, result
(0, 0x02)
```
"""
function project!(stabilizer::Stabilizer,pauli::PauliOperator;keep_result::Bool=true,phases::Bool=true)
    anticommutes = 0
    n = size(stabilizer,1)
    for i in 1:n
        if comm(pauli,stabilizer,i)!=0x0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        if keep_result
            canonicalize!(stabilizer; phases=phases)
            gen = generate!(copy(pauli), stabilizer, phases=phases)
            result = isnothing(gen) ? nothing : gen[1].phase[]
        else
            result = nothing
        end
    else
        for i in anticommutes+1:n
            if comm(pauli,stabilizer,i)!=0
                rowmul!(stabilizer, i, anticommutes; phases=phases)
            end
        end
        stabilizer[anticommutes] = pauli
        result = nothing
    end
    stabilizer, anticommutes, result
end

function project!(d::Destabilizer,pauli::PauliOperator;keep_result::Bool=true,phases::Bool=true)
    anticommutes = 0
    stabilizer = stabilizerview(d)
    destabilizer = destabilizerview(d)
    n = size(stabilizer,1)
    for i in 1:n
        if comm(pauli,stabilizer,i)!=0x0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        if n != nqubits(stabilizer)
            throw(BadDataStructure("`Destabilizer` can not efficiently (faster than n^3) detect whether you are projecting on a stabilized or a logical operator. Switch to one of the `Mixed*` data structures.",
                                   :project!,
                                   :Destabilizer))
        end
        if keep_result
            new_pauli = zero(pauli)
            for i in 1:n
                if comm(pauli,destabilizer,i)!=0
                    # TODO, this is just a long explicit way to write it... learn more about broadcast
                    phases && (new_pauli.phase[] = prodphase(stabilizer, new_pauli, i))
                    new_pauli.xz .⊻= @view stabilizer.xzs[i,:]
                end
            end
            result = new_pauli.phase[]
        else
            result = nothing
        end
    else
        for i in anticommutes+1:n
            if comm(pauli,stabilizer,i)!=0
                rowmul!(stabilizer, i, anticommutes; phases=phases)
            end
        end
        for i in 1:n
            if i!=anticommutes && comm(pauli,destabilizer,i)!=0
                rowmul!(d.tab, i, n+anticommutes; phases=false)
            end
        end
        destabilizer[anticommutes] = stabilizer[anticommutes]
        stabilizer[anticommutes] = pauli
        result = nothing
    end
    d, anticommutes, result
end

function project!(ms::MixedStabilizer,pauli::PauliOperator;keep_result::Bool=true,phases::Bool=true)
    _, anticom_index, res = project!(stabilizerview(ms), pauli; keep_result=keep_result, phases=phases)
    if anticom_index==0 && isnothing(res)
        ms.tab[ms.rank+1] = pauli
        if keep_result
            ms.rank += 1
        else
            canonicalize!(ms.tab[1:ms.rank+1]; phases=phases)
            if ~all(==(UInt64(0)), @view ms.tab.xzs[ms.rank+1,:])
                ms.rank += 1
            end
        end
    end
    ms, anticom_index, res # TODO CHECK THIS res
end

function anticomm_update_rows(tab,pauli,r,n,anticommutes,phases) # TODO Ensure there are no redundant `comm` checks that can be skipped
    chunks = size(tab.xzs,2)
    for i in r+1:n
        if comm(pauli,tab,i)!=0
            rowmul!(tab, i, n+anticommutes; phases=phases)
        end
    end
    for i in n+anticommutes+1:2n
        if comm(pauli,tab,i)!=0
            rowmul!(tab, i, n+anticommutes; phases=phases)
        end
    end
    for i in 1:r
        if i!=anticommutes && comm(pauli,tab,i)!=0
            rowmul!(tab, i, n+anticommutes; phases=false)
        end
    end
end

function project!(d::MixedDestabilizer,pauli::PauliOperator;keep_result::Bool=true,phases::Bool=true)
    anticommutes = 0
    tab = d.tab
    stabilizer = stabilizerview(d)
    destabilizer = destabilizerview(d)
    r = d.rank
    n = nqubits(d)
    for i in 1:r # TODO use something like findfirst
        if comm(pauli,stabilizer,i)!=0x0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        anticomlog = 0
        for i in r+1:n # TODO use something like findfirst
            if comm(pauli,tab,i)!=0x0
                anticomlog = i
                break
            end
        end
        if anticomlog==0
            for i in n+r+1:2*n # TODO use something like findfirst
                if comm(pauli,tab,i)!=0x0
                    anticomlog = i
                    break
                end
            end
        end
        if anticomlog!=0
            if anticomlog<=n
                rowswap!(tab, r+1+n, anticomlog)
                n!=r+1 && anticomlog!=r+1 && rowswap!(tab, r+1, anticomlog+n)
            else
                rowswap!(tab, r+1, anticomlog-n)
                rowswap!(tab, r+1+n, anticomlog)
            end
            anticomm_update_rows(tab,pauli,r+1,n,r+1,phases)
            d.rank += 1
            tab[r+1] = tab[n+r+1]
            tab[n+r+1] = pauli
            result = nothing
        else
            if keep_result
                new_pauli = zero(pauli)
                for i in 1:r
                    if comm(pauli,destabilizer,i)!=0
                        # TODO, this is just a long explicit way to write it... learn more about broadcast
                        phases && (new_pauli.phase[] = prodphase(stabilizer, new_pauli, i))
                        new_pauli.xz .⊻= @view stabilizer.xzs[i,:]
                    end
                end
                result = new_pauli.phase[]
            else
                result = nothing
            end
        end
    else
        anticomm_update_rows(tab,pauli,r,n,anticommutes,phases)
        destabilizer[anticommutes] = stabilizer[anticommutes]
        stabilizer[anticommutes] = pauli
        result = nothing
    end
    d, anticommutes, result
end

"""
Trace out a qubit.
""" # TODO all of these should raise an error if length(qubits)>rank
function traceout!(s::MixedStabilizer, qubits::AbstractVector{T}; phases=true) where {T<:Integer} # TODO implement it on the other state data structures.
    _,i = canonicalize_rref!(s,qubits;phases=phases)
    s.rank = i
end

function traceout!(s::MixedDestabilizer, qubits::AbstractVector{T}; phases=true) where {T<:Integer} # TODO implement it on the other state data structures.
    _,i = canonicalize_rref!(s,qubits;phases=phases)
    s.rank = i
end

##############################
# Unitary Clifford Operations
##############################

function Base.:(*)(p::AbstractCliffordOperator, s::AbstractStabilizer; phases::Bool=true)
    s = copy(s)
    apply!(s,p; phases=phases)
end

function apply!(s::Stabilizer, p::PauliOperator; phases::Bool=true)
    phases || return s
    for i in eachindex(s)
        s.phases[i] = (s.phases[i]+comm(p,s,i)<<1+p.phase[]<<1)&0x3
    end
    s
end

"""
Clifford Operator specified by the mapping of the basis generators.

```jldoctest
julia> CNOT
X_ ⟼ + XX
_X ⟼ + _X
Z_ ⟼ + Z_
_Z ⟼ + ZZ

julia> phase_gate = C"Y
                      Z"
X ⟼ + Y
Z ⟼ + Z

julia> stab = S"XI
                IZ";

julia> entangled = CNOT*stab
+ XX
+ ZZ
```
"""
struct CliffordOperator{Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{UInt64}} <: AbstractCliffordOperator
    phases::Tv
    nqubits::Int
    xztox::Tm
    xztoz::Tm
end

function CliffordOperator(paulis::AbstractVector{PauliOperator{Tz,Tv}}) where {Tz<:AbstractArray{UInt8,0},Tv<:AbstractVector{UInt64}}
    xztox = vcat((p.xbit' for p in paulis)...)'
    xztoz = vcat((p.zbit' for p in paulis)...)'
    xztox = vcat((vcat(BitArray(xztox[i,1:end÷2]).chunks,BitArray(xztox[i,end÷2+1:end]).chunks)'
                 for i in 1:size(xztox,1))...)
    xztoz = vcat((vcat(BitArray(xztoz[i,1:end÷2]).chunks,BitArray(xztoz[i,end÷2+1:end]).chunks)'
                 for i in 1:size(xztoz,1))...)
    CliffordOperator(vcat((p.phase for p in paulis)...), paulis[1].nqubits, xztox, xztoz)
end

macro C_str(a)
    paulis = [eval(quote @P_str($(strip(s))) end) for s in split(a,'\n')] #TODO seriously!?
    CliffordOperator(paulis)
end

Base.:(==)(l::CliffordOperator, r::CliffordOperator) = r.nqubits==l.nqubits && r.phases==l.phases && r.xztox==l.xztox && r.xztoz==l.xztoz

function clifford_transpose(c::CliffordOperator)
    n = c.nqubits
    xtoxs = []
    ztozs = []
    ztoxs = []
    xtozs = []
    for r in 1:n
        xtox = BitArray(UndefInitializer(),(n,))
        ztoz = BitArray(UndefInitializer(),(n,))
        xtoz = BitArray(UndefInitializer(),(n,))
        ztox = BitArray(UndefInitializer(),(n,))
        xtox.chunks = c.xztox[r,1:end÷2]
        ztoz.chunks = c.xztoz[r,end÷2+1:end]
        xtoz.chunks = c.xztoz[r,1:end÷2]
        ztox.chunks = c.xztox[r,end÷2+1:end]
        push!(xtoxs,xtox')
        push!(ztozs,ztoz')
        push!(xtozs,xtoz')
        push!(ztoxs,ztox')
    end
    xtoxs = vcat(xtoxs...)
    ztozs = vcat(ztozs...)
    xtozs = vcat(xtozs...)
    ztoxs = vcat(ztoxs...)
    xtoxs, ztozs, xtozs, ztoxs
end

function getallpaulis_(c::CliffordOperator)
    xtoxs, ztozs, xtozs, ztoxs = clifford_transpose(c)
    ops = PauliOperator{Array{UInt8,0},Vector{UInt64}}[] # TODO is it really necessary to specify the type this precisely!?
    for i in 1:2*c.nqubits
        if i>c.nqubits
            push!(ops,PauliOperator(c.phases[i], ztoxs[:,i-c.nqubits], ztozs[:,i-c.nqubits]))
        else
            push!(ops,PauliOperator(c.phases[i], xtoxs[:,i], xtozs[:,i]))
        end
    end
    ops
end

function Base.getindex(c::CliffordOperator, i::Int)
    xtoxs, ztozs, xtozs, ztoxs = clifford_transpose(c)
    if i>c.nqubits
        PauliOperator(c.phases[i], ztoxs[:,i-c.nqubits], ztozs[:,i-c.nqubits])
    else
        PauliOperator(c.phases[i], xtoxs[:,i], xtozs[:,i])
    end
end

function Base.show(io::IO, c::CliffordOperator)
    xtoxs, ztozs, xtozs, ztoxs = clifford_transpose(c)
    n = c.nqubits
    for i in 1:n
        print(io, repeat("_",i-1),"X",repeat("_",n-i), " ⟼ ")
        print(io, ["+ ","+i","- ","-i"][c.phases[i]+1])
        print(io, xz2str(xtoxs[:,i],xtozs[:,i]))
        println(io)
    end
    for i in 1:n
        print(io, repeat("_",i-1),"Z",repeat("_",n-i), " ⟼ ")
        print(io, ["+ ","+i","- ","-i"][c.phases[i+n]+1])
        print(io, xz2str(ztoxs[:,i],ztozs[:,i]))
        println(io)
    end
end

function Base.copy(c::CliffordOperator)
    CliffordOperator(copy(c.phases),c.nqubits,copy(c.xztox),copy(c.xztoz))
end

# TODO create Base.permute! and getindex(..., permutation_array)
function permute(c::CliffordOperator,p::AbstractArray{T,1} where T) # TODO this is extremely slow stupid implementation
    ops = getallpaulis_(c)
    CliffordOperator([ops[i][p] for i in 1:2*c.nqubits][vcat(p,p.+c.nqubits)])
end

function apply!(s::Stabilizer, c::CliffordOperator; phases::Bool=true)
    new_stabrowx = zero(s.xzs[1,1:end÷2])
    new_stabrowz = zero(s.xzs[1,1:end÷2])
    xztox = zero(c.xztox[1,:])
    xztoz = zero(c.xztox[1,:])
    phase_flips = zero(xztoz[1:end÷2])
    for row_stab in eachindex(s)
        fill!(new_stabrowx, zero(eltype(new_stabrowx)))
        fill!(new_stabrowz, zero(eltype(new_stabrowz)))
        for row_clif in 1:nqubits(s)
            bigrow = _div64(row_clif-1)+1
            smallrow = _mod64(row_clif-1)
            @inbounds @simd for i in 1:length(xztox)
                xztox[i] = c.xztox[row_clif,i] & s.xzs[row_stab,i]
                xztoz[i] = c.xztoz[row_clif,i] & s.xzs[row_stab,i]
            end
            new_stabrowx[bigrow] |= xor_bits_(reduce(⊻,xztox)) << smallrow
            new_stabrowz[bigrow] |= xor_bits_(reduce(⊻,xztoz)) << smallrow
            phases && (phase_flips .= (@view xztoz[1:end÷2]) .& (@view xztox[end÷2+1:end]))
            phases && (s.phases[row_stab] = (s.phases[row_stab]+sum(count_zeros, phase_flips)<<1)&0x3)
        end
        s.xzs[row_stab,1:end÷2] = new_stabrowx
        s.xzs[row_stab,end÷2+1:end] = new_stabrowz
    end
    s
end

function apply!(s::Stabilizer, c::CliffordOperator, single_qbit_offset::Int)
    bigs = _div64(single_qbit_offset-1)+1
    smalls = _mod64(single_qbit_offset-1)
    lowbit = UInt64(0x1)
    for row_stab in eachindex(s)
        xztox = (c.xztox[1,1] & (s.xzs[row_stab,bigs]>>smalls)) ⊻ (c.xztox[1,2] & (s.xzs[row_stab,bigs+end÷2]>>smalls))
        xztoz = (c.xztoz[1,1] & (s.xzs[row_stab,bigs]>>smalls)) ⊻ (c.xztoz[1,2] & (s.xzs[row_stab,bigs+end÷2]>>smalls))
        s.phases[row_stab] = (s.phases[row_stab]+count_zeros(xztoz & xztox)<<1)&0x3

        s.xzs[row_stab,bigs] &= ~(lowbit<<smalls)
        s.xzs[row_stab,bigs] |= (xztox<<smalls)
        s.xzs[row_stab,end÷2+bigs] &= ~(lowbit<<smalls)
        s.xzs[row_stab,end÷2+bigs] |= (xztoz<<smalls)
    end
    s
end

function tensor_pow(op::CliffordOperator,power::Integer,mem::Dict{Integer,CliffordOperator})
    if power==1
        return op
    elseif haskey(mem,power)
        return mem[power]
    end
    half,rest = divrem(power,2)
    phalf = get!(mem,half) do
        tensor_pow(op,half,mem)
    end
    res = phalf⊗phalf
    if rest!=0
        prest = get!(mem,rest) do
            tensor_pow(op,half,mem)
        end
        res = res⊗prest
    end
    res
end

function tensor_pow(op::CliffordOperator,power::Integer)
    tensor_pow(op,power,Dict{Integer,CliffordOperator}())
end

const CNOT = C"XX
               IX
               ZI
               ZZ"

const SWAP = C"IX
               XI
               IZ
               ZI"

const Hadamard = C"Z
                   X"

const Phase = C"Y
                Z"

const CliffordId = C"X
                     Z"

##############################
# Helpers for Clifford Operators
##############################

function (⊗)(l::CliffordOperator, r::CliffordOperator) # TODO this is extremely slow stupid implementation
    opsl = getallpaulis_(l)
    opsr = getallpaulis_(r)
    onel = zero(opsl[1])
    oner = zero(opsr[1])
    opsl = [l⊗oner for l in opsl]
    opsr = [onel⊗r for r in opsr]
    CliffordOperator(vcat(opsl[1:end÷2],opsr[1:end÷2],opsl[end÷2+1:end],opsr[end÷2+1:end]))
end

function Base.:(*)(l::AbstractCliffordOperator, r::CliffordOperator)
    rstab = Stabilizer(getallpaulis_(r)) # TODO this is a bit awkward... and fragile... turning a CliffordOp into a Stabilizer
    apply!(rstab,l)
    CliffordOperator([rstab[i] for i in eachindex(rstab)])
end

##############################
# Helpers for binary codes
##############################

function stab_to_gf2(s::Stabilizer)
    xbits = vcat((s[i].xbit' for i in eachindex(s))...)
    zbits = vcat((s[i].zbit' for i in eachindex(s))...)
    H = hcat(xbits,zbits)
end

function gf2_gausselim!(H) # equivalent to just taking the canonicalized stabilizer
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

function gf2_isinvertible(H) # TODO can be smarter and exit earlier. And should check for squareness.
    ut = gf2_gausselim!(copy(H))
    all((ut[i, i] for i in 1:size(ut,1)))
end

function gf2_invert(H)
    id = zero(H)
    s = size(H,1)
    for i in 1:s id[i,i]=true end
    M = hcat(H,id)
    gf2_gausselim!(M)
    M[:,s+1:end]
end

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
    G[:,perm_inverse(sindx)]
end

function perm_inverse(perm)
    [findfirst(a->l==a,perm) for l in 1:length(perm)]
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
# Common objects
##############################

function single_z(n,i)
    xs = falses(n)
    zs = falses(n)
    zs[i] = true
    PauliOperator(0x0,xs,zs)
end

function single_x(n,i)
    xs = falses(n)
    zs = falses(n)
    xs[i] = true
    PauliOperator(0x0,xs,zs)
end

Base.zero(::Type{PauliOperator}, n) = PauliOperator(0x0,falses(n),falses(n))
Base.zero(p::PauliOperator) = PauliOperator(0x0,falses(p.nqubits),falses(p.nqubits))
Base.zero(::Type{Stabilizer}, n, m) = Stabilizer(zeros(UInt8,n),falses(n,m),falses(n,m))
Base.zero(::Type{Stabilizer}, n) = Stabilizer(zeros(UInt8,n),falses(n,n),falses(n,n))
Base.zero(s::Stabilizer) = Stabilizer(zeros(UInt8,size(s,1)),falses(size(s)...),falses(size(s)...))

function Base.one(::Type{Stabilizer}, n; basis=:Z)
    if basis==:X
        Stabilizer(LinearAlgebra.I(n),falses(n,n))
    elseif basis==:Y
        Stabilizer(LinearAlgebra.I(n),LinearAlgebra.I(n))
    elseif basis==:Z
        Stabilizer(falses(n,n),LinearAlgebra.I(n))
    else
        throw(ErrorException("`basis` should be one of :X, :Y, or :Z"))
    end
end
Base.one(s::Stabilizer; basis=:Z) = one(Stabilizer, nqubits(s); basis=basis)
Base.one(::Type{Destabilizer}, n) = Destabilizer(vcat(one(Stabilizer, n, basis=:X),one(Stabilizer, n, basis=:Z)),noprocessing=true)
function Base.one(::Type{MixedDestabilizer}, r, n)
    d = one(Stabilizer, n; basis=:X)
    s = one(Stabilizer, n; basis=:Z)
    MixedDestabilizer(vcat(d,s),r)
end

##############################
# Random objects
##############################

random_pauli(n; nophase=false) = PauliOperator(nophase ? 0x0 : rand(0x0:0x3), rand(Bool,n), rand(Bool,n))
function random_pauli(n,p; nophase=false,nonidbranch=false)
    x = falses(n)
    z = falses(n)
    if nonidbranch
        definite = rand(1:n)
        p=(p/(1-(1-2p)^n) - 1/n/3)*n/(n-1)
    end
    for i in 1:n
        r = rand()
        if nonidbranch && definite==i
            r *= 3p
        end
        if (r<=2p) x[i]=true end
        if (p<r<=3p) z[i]=true end
    end
    PauliOperator(nophase ? 0x0 : rand(0x0:0x3), x, z)
end

function random_invertible_gf2(n)
    while true
        mat = rand(Bool,n,n)
        gf2_isinvertible(mat) && return mat
    end
end

# function random_cnot_clifford(n) = ... #TODO

function random_stabilizer(n) # TODO this is vaguelly based on an unsupported slide deck off the internet. Probably incorrectly implemented too.
    cx = falses(n,n)
    cz = falses(n,n)
    for i in 1:n
        cx[i,i], cz[i,i] = rand([(true,true),(true,false),(false,true)])
    end
    C = random_invertible_gf2(n)
    CinvT = gf2_invert(C)'
    cx = Bool.((cx * C) .% 2)
    cz = Bool.((cz * CinvT) .% 2)
    Stabilizer(rand([0x0,0x2],n), cx, cz)
end

random_stabilizer(r,n) = random_stabilizer(n)[Random.randperm(n)[1:r]]

function random_singlequbitop(n)
    xtox = [falses(n) for i in 1:n]
    ztox = [falses(n) for i in 1:n]
    xtoz = [falses(n) for i in 1:n]
    ztoz = [falses(n) for i in 1:n]
    for i in 1:n
        gate = rand(1:6)
        if gate<5
            xtox[i][i] = true
            xtoz[i][i] = true
            ztox[i][i] = true
            ztoz[i][i] = true
            [xtox,ztox,xtoz,ztoz][gate][i][i] = false
        elseif gate==5
            xtox[i][i] = true
            ztoz[i][i] = true
        else
            xtoz[i][i] = true
            ztox[i][i] = true
        end
    end
    c = CliffordOperator(zeros(UInt8,n*2), n,
                         vcat((vcat(x2x.chunks,z2x.chunks)' for (x2x,z2x) in zip(xtox,ztox))...),
                         vcat((vcat(x2z.chunks,z2z.chunks)' for (x2z,z2z) in zip(xtoz,ztoz))...)
        )
end

##############################
# Plotting
##############################

RecipesBase.@recipe function f(s::Stabilizer; xzcomponents=:split)
    seriestype  := :heatmap
    aspectratio := :equal
    yflip := true
    colorbar := false
    grid := false
    framestyle := :none
    if xzcomponents==:split
        stab_to_gf2(s)
    elseif xzcomponents==:together
        h = stab_to_gf2(s)
        h[:,1:end÷2]*2 + h[:,end÷2+1:end]
    else
        throw(ErrorException("`xzcomponents` should be `:split` or `:together`"))
    end
end

end #module
