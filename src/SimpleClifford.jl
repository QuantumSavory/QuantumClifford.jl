"""
A module for simulation of Clifford circuits.
It does not yet use the "destabilizer" tableau formalism, hence many algorithms are
cubic instead of quadratic.

Single-qubit Paulis will be encoded as elements of F(2,2) in the low two bits of an `UInt8`.
"""
module SimpleClifford

# TODO remove the triple quotes

# TODO is this https://www.cs.umd.edu/~amchilds/teaching/w10/project-sample.pdf useful?

# TODO the traceout and reset functions can be significantly optimized (they require a ridiculous repetition of canonicalizations and projections right now).

# TODO do not mix up getindex and view (currently getindex is sometimes a view and there is no official view)

# Operations between Clifford operators are very slow

export @P_str, PauliOperator, ⊗, I, X, Y, Z, permute,
    @S_str, Stabilizer, prodphase, comm, ⊕, isfullstabilizer, canonicalize!,
    generate!, project!, reset_qubits!, traceout_qubits!,
    apply!,
    CliffordOperator, @C_str, CNOT, SWAP, Hadamard, Phase, CliffordId

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

julia> p.phase, p.xbit, p.zbit
(0x02, Bool[0, 0, 1, 1], Bool[0, 1, 0, 1])

julia> p.xz
2-element Array{UInt64,1}:
 0x000000000000000c
 0x000000000000000a
```
"""
struct PauliOperator{Tz<:AbstractArray{UInt8,0},Tv<:AbstractVector{UInt64}} <: AbstractCliffordOperator
    phase::Tz
    nqbits::Int
    xz::Tv
end

PauliOperator(phase::UInt8, nqbits::Int, xz::Tv) where Tv<:AbstractVector{UInt64} = PauliOperator(fill(phase,()), nqbits, xz)
PauliOperator(phase::UInt8, x::T, z::T) where T<:AbstractVector{Bool} = PauliOperator(fill(phase,()), length(x), vcat(BitVector(x).chunks,BitVector(z).chunks))

function Base.getproperty(p::PauliOperator, name::Symbol)
    if name==:xview
        @view p.xz[1:end÷2]
    elseif name==:zview
        @view p.xz[end÷2+1:end]
    elseif name==:xbit
        b = BitArray(UndefInitializer(),(p.nqbits,))
        b.chunks = p.xview
        b
    elseif name==:zbit
        b = BitArray(UndefInitializer(),(p.nqbits,))
        b.chunks = p.zview
        b
    else
        getfield(p, name)
    end
end

Base.propertynames(p::PauliOperator, private=false) = (:phases,:nqbits,:xz,:xbit,:zbit,:xview,:zview)

macro P_str(a)
    letters = filter(x->occursin(x,"_IZXY"),a)
    phase = phasedict[strip(filter(x->!occursin(x,"_IZXY"),a))]
    PauliOperator(phase, [l=='X'||l=='Y' for l in letters], [l=='Z'||l=='Y' for l in letters])
end

Base.size(pauli::PauliOperator) = pauli.nqbits

xz2str(x,z) = join(toletter[e] for e in zip(x,z))

Base.show(io::IO, p::PauliOperator) = print(io, ["+ ","+i","- ","-i"][p.phase[]+1]*xz2str(p.xbit,p.zbit))

Base.:(==)(l::PauliOperator, r::PauliOperator) = r.phase==l.phase && r.nqbits==l.nqbits && r.xz==l.xz

Base.hash(p::PauliOperator, h::UInt) = hash((p.phase,p.nqbits,p.xz), h)

Base.copy(p::PauliOperator) = PauliOperator(copy(p.phase),p.nqbits,copy(p.xz))

Base.one(p::PauliOperator) = PauliOperator(zeros(UInt8),p.nqbits,zero(p.xz))

##############################
# Pauli Operator Helpers
##############################

@inline function prodphase_(x1::AbstractVector{UInt64},z1::AbstractVector{UInt64},x2::AbstractVector{UInt64},z2::AbstractVector{UInt64})::UInt64
    pos = (.~z2 .& x2 .& .~x1 .& z1) .| (z2 .& .~x2 .& x1 .& z1) .| (z2 .&   x2 .& x1 .& .~z1)
    neg = (  z2 .& x2 .& .~x1 .& z1) .| (.~z2 .& x2 .& x1 .& z1) .| (z2 .& .~x2 .& x1 .& .~z1)
    unsigned(sum(count_ones,pos) - sum(count_ones,neg))#<<1
end

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
@inline prodphase(l::PauliOperator, r::PauliOperator)::UInt8 = (l.phase+r.phase+prodphase_(l.xview,l.zview,r.xview,r.zview))&0x3

@inline function xor_bits_(v::UInt64)
    v ⊻= v >> 32
    v ⊻= v >> 16
    v ⊻= v >> 8
    v ⊻= v >> 4
    v ⊻= v >> 2
    v ⊻= v >> 1
    return v&1
end

@inline function comm_(
        x1::AbstractVector{UInt64},
        z1::AbstractVector{UInt64},
        x2::AbstractVector{UInt64},
        z2::AbstractVector{UInt64})::UInt8 # based on the twisted product from arxiv 0304161 or arxiv 9608006
    xor_bits_(reduce(⊻,((z1 .& x2) .⊻ (z2 .& x1))))
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
@inline comm(l::PauliOperator, r::PauliOperator)::UInt8 = comm_(l.xz[1:end÷2],l.xz[end÷2+1:end],r.xz[1:end÷2],r.xz[end÷2+1:end])
#comm(l::PauliOperator, r::PauliOperator)::UInt8 = comm_(l.xview,l.zview,r.xview,r.zview)


function Base.:(*)(l::PauliOperator, r::PauliOperator)
    PauliOperator(prodphase(l,r), l.nqbits, l.xz .⊻ r.xz)
end

(⊗)(l::PauliOperator, r::PauliOperator) = PauliOperator((l.phase+r.phase)&0x3, vcat(l.xbit,r.xbit) , vcat(l.zbit,r.zbit))

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

# TODO create Base.permute! and getindex(..., permutation_array)
function permute(p::PauliOperator,perm::AbstractArray{T,1} where T)
    PauliOperator(p.phase[],p.xbit[perm],p.zbit[perm])
end

const I = P"I"
const Z = P"Z"
const X = P"X"
const Y = P"Y"

##############################
# Stabilizers
##############################

"""
Stabilizer, i.e. a list of commuting multi-qubit Hermitian Pauli operators.

Instances can be created with the `S` custom string macro or
as direct sum of other stabilizers.

```jldoctest stabilizer
julia> s = S\"""XXX
               ZZI
               IZZ\"""
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
commutativity of all rows, hermiticity of all rows).

See also: [`PauliOperator`](@ref), [`canonicalize!`](@ref)
"""
struct Stabilizer{Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{UInt64}}
    phases::Tv
    nqbits::Int
    xzs::Tm
end

Stabilizer(paulis::AbstractVector{PauliOperator{Tz,Tv}}) where {Tz<:AbstractArray{UInt8,0},Tv<:AbstractVector{UInt64}} = Stabilizer(vcat((p.phase for p in paulis)...), paulis[1].nqbits, vcat((p.xz' for p in paulis)...))

macro S_str(a)
    paulis = [eval(quote @P_str($(strip(s))) end) for s in split(a,'\n')] #TODO seriously!?
    Stabilizer(paulis)
end

Base.getindex(stab::Stabilizer, i::Int) = PauliOperator((@view stab.phases[i]), stab.nqbits, (@view stab.xzs[i,:]))
Base.getindex(stab::Stabilizer, r::UnitRange) = Stabilizer((@view stab.phases[r]), stab.nqbits, (@view stab.xzs[r,:]))

function Base.setindex!(stab::Stabilizer, pauli::PauliOperator, i)
    stab.phases[i] = pauli.phase[]
    stab.xzs[i,:] = pauli.xz
    pauli
end

Base.firstindex(stab::Stabilizer) = 1

Base.lastindex(stab::Stabilizer) = length(stab.phases)

Base.show(io::IO, s::Stabilizer) = print(io,
                                         join([["+ ","+i","- ","-i"][s[i].phase[]+1]*xz2str(s[i].xbit,s[i].zbit)
                                               for i in firstindex(s):lastindex(s)],
                                              '\n'))

Base.:(==)(l::Stabilizer, r::Stabilizer) = r.nqbits==l.nqbits && r.phases==l.phases && r.xzs==l.xzs

Base.hash(s::Stabilizer, h::UInt) = hash(s.nqbits, s.phases, s.xzs, h)

Base.copy(s::Stabilizer) = Stabilizer(copy(s.phases), s.nqbits, copy(s.xzs))

@inline function rowswap!(s::Stabilizer, i, j; phases::Bool=true) # Written only so we can avoid copying in `canonicalize!`
    (i == j) && return
    phases && begin s.phases[i], s.phases[j] = s.phases[j], s.phases[i] end
    @simd for k in 1:size(s.xzs,2)
        s.xzs[i,k], s.xzs[j,k] = s.xzs[j,k], s.xzs[i,k]
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

Assumes all operators do commute, but it permits redundant generators or
incomplete list of generators in the input.

It permits non-Hermitian generators, even though it is physically meaningless.

Based on arxiv:1210.6646.

```jldoctest
julia> ghz = S\"""XXXX
                 ZZII
                 IZZI
                 IIZZ\""";

julia> canonicalize!(ghz)
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ

julia> canonicalize!(S\"""XXXX
                         IZZI
                         IIZZ\""")
+ XXXX
+ _Z_Z
+ __ZZ

julia> canonicalize!(S\"""XXXX
                         ZZII
                         IZZI
                         IZIZ
                         IIZZ\""")
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ
+ ____
```

See arxiv:0505036 for other types of canonicalization.
"""
function canonicalize!(stabilizer::Stabilizer; phases::Bool=true) # TODO simplify by using the new Pauli like interface instead of f22array
    xzs = stabilizer.xzs
    xs = @view xzs[:,1:end÷2]
    zs = @view xzs[:,end÷2+1:end]
    lowbit = UInt64(0x1)
    zero64 = UInt64(0x0)
    rows = length(stabilizer.phases)
    columns = stabilizer.nqbits
    i = 1
    for j in 1:columns
        # find first row with X or Y in col `j`
        jbig = j÷64+1  # TODO use _div and _mod
        jsmall = lowbit<<((j-1)%64)  # TODO use _div and _mod
        k = findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                      xs[i:end,jbig])
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases=phases)
            for m in 1:rows
                if xs[m,jbig]&jsmall!=zero64 && m!=i # if X or Y
                    xzs[m,:] .⊻= xzs[i,:]
                    phases && (stabilizer.phases[m] = (stabilizer.phases[m]+stabilizer.phases[i]+prodphase_(xs[m,:],zs[m,:],xs[i,:],zs[i,:]))&0x3)
                end
            end
            i += 1
        end
    end
    for j in 1:columns
        # find first row with Z in col `j`
        jbig = j÷64+1  # TODO use _div and _mod
        jsmall = lowbit<<((j-1)%64)  # TODO use _div and _mod
        k = findfirst(e->e&(jsmall)!=zero64,
                      zs[i:end,jbig])
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases=phases)
            for m in 1:rows
                if zs[m,jbig]&jsmall!=zero64 && m!=i # if Z or Y
                    xzs[m,:] .⊻= xzs[i,:]
                    phases && (stabilizer.phases[m] = (stabilizer.phases[m]+stabilizer.phases[i]+prodphase_(xs[m,:],zs[m,:],xs[i,:],zs[i,:]))&0x3)
                end
            end
            i += 1
        end
    end
    stabilizer
end

function ishermitian() # TODO write it both for paulis and stabilizers (ugh... stabilizers are always so)
end

function isfullstabilizer(stabilizer::Stabilizer) # TODO update
#=    s = stabilizer.F22array
    n,m = size(s)
    if n!=m
        return false
    end
    for i in 1:n
        for j in i+1:n
            if comm(s[i,:],s[j,:])!=0x0
                return false
            end
        end
    end
    return true=#
end

function ⊕(l::Stabilizer, r::Stabilizer)
    lone = one(l[1])
    rone = one(r[1])
    paulis = vcat([l[i]⊗rone for i in firstindex(l):lastindex(l)],
                  [lone⊗r[i] for i in firstindex(r):lastindex(r)]
                 )
    Stabilizer(paulis)
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
julia> ghz = S\"""XXXX
                 ZZII
                 IZZI
                 IIZZ\""";

julia> canonicalize!(ghz)
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ

julia> generate!(P"-ZIZI", ghz)
(- ____, [2, 4])
```
"""
function generate!(pauli::PauliOperator, stabilizer::Stabilizer; phases::Bool=true) # TODO there is stuff that can be abstracted away here and in canonicalize!
    rows = length(stabilizer.phases)
    columns = stabilizer.nqbits
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
        jbig = i÷64+1  # TODO use _div and _mod
        jsmall = lowbit<<((i-1)%64)  # TODO use _div and _mod
        used += findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                          xs[used+1:end,jbig])
        # TODO, this is just a long explicit way to write it... learn more about broadcast
        phases && (pauli.phase[] = prodphase(pauli,stabilizer[used]))
        pauli.xz .⊻= xzs[used,:]
        push!(used_indices, used)
    end
    # remove Zs
    while (i=unsafe_bitfindnext_(pz,1)) !== nothing
        jbig = i÷64+1  # TODO use _div and _mod
        jsmall = lowbit<<((i-1)%64)  # TODO use _div and _mod
        used += findfirst(e->e&jsmall!=zero64, # TODO some form of reinterpret might be faster than equality check
                          zs[used+1:end,jbig])
        # TODO, this is just a long explicit way to write it... learn more about broadcast
        phases && (pauli.phase[] = prodphase(pauli,stabilizer[used]))
        pauli.xz .⊻= xzs[used,:]
        push!(used_indices, used)
    end
    pauli, used_indices
end

"""
Project the state of a Stabilizer on the two eigenspaces of a Pauli operator.

Assumes non-identity operators with phase ±1 (they need to be Hermitian).
The projection is done inplace.

It returns a 

 - non-canonicalized stabilizer
 - the index of the row where the non-commuting operator was (that row is now equal to `pauli` and has an uncertain phase)
 - and the result of the projection if there was no non-cummuting operator

If `keep_result==false` that result of the projection is not computed,
sparing a canonicalization operation.

Here is an example of a projection destroing entanglement:

```jldoctest
julia> ghz = S\"""XXXX
                 ZZII
                 IZZI
                 IIZZ\""";

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
julia> s = S\"""ZII
               IXI
               IIY\""";

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

The measurement result can be `nothing` or the `anticom_index` can be `0`.
This will turn into a neater interface at some point.
"""# TODO make a neater interface as suggested just above.
function project!(stabilizer::Stabilizer,pauli::PauliOperator;keep_result::Bool=true,phases::Bool=true)
    # @assert !all(pauli.F22array .== I)
    # @assert pauli.phase ∈ [0x0, 0x2]
    anticommutes = 0                                           
    n = length(stabilizer.phases)
    for i in 1:n
        if comm(pauli,stabilizer[i])!=0x0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        if keep_result
            canonicalize!(stabilizer; phases=phases)
            new_pauli, _ = generate!(copy(pauli), stabilizer, phases=phases)
            result = new_pauli.phase
        else
            result = nothing
        end
    else
        for i in anticommutes+1:n
            if comm(pauli,stabilizer[i])!=0
                # TODO, this is just a long explicit way to write it... learn more about broadcast
                phases && (stabilizer.phases[i] = prodphase(stabilizer[i], stabilizer[anticommutes]))
                stabilizer.xzs[i,:] .⊻= stabilizer.xzs[anticommutes,:]
            end
        end
        stabilizer[anticommutes] = pauli
        result = nothing
    end
    stabilizer, anticommutes, result
end

#="""
Reset specified qubits to have a new stabilizer state.

The qubits are assumed to be unentangled
and the stabilizer is assumed canonicalized.

```jldoctest
julia> s_entangled = S"XXI
                       ZZI
                       IIZ";

julia> s_product = S"ZII
                     IXI
                     IIY";

julia> newstab = S"-YX
                   +XY";

julia> reset_qubits!(s_product, newstab, [1,3])
- Y_X
+ _X_
+ X_Y

julia> reset_qubits!(s_entangled, newstab, [1,3])
ERROR: AssertionError: the qubits to be reset are entangled
[...]
```
"""=# # TODO fix/update this
function reset_qubits!(stabilizer::Stabilizer, newsubstabilizer::Stabilizer, qubits) # TODO type the qubits as an index
#=    s = stabilizer.F22array
    origrows, origcols = size(stabilizer)
    rows = zeros(Bool, origrows)
    for q in enumerate(qubits)
        rows .|= _I.!=s[:,q]
    end
    @assert sum(rows) == length(qubits) "the qubits to be reset are entangled"
    # TODO the check above is not enough
    stabilizer.phasesF22array[rows,[1,[q+1 for q in qubits]...]] = newsubstabilizer.phasesF22array # TODO nicer api
    stabilizer =#
end

#="""
Trace out qubits.

The qubits are assumed to be unentangled
and the stabilizer is assumed canonicalized.

```jldoctest
julia> s_entangled = S"XXI
                       ZZI
                       IIZ";

julia> s_product = S"ZII
                     IXI
                     IIY";

julia> traceout_qubits!(s_product, [1,3])
+ X

julia> traceout_qubits!(s_entangled, [1,3])
ERROR: AssertionError: the qubits to be reset are entangled
[...]
```
"""=# # TODO fix/update this # TODO it is not exactly mutable, rather it mutates and then returns a view...
function traceout_qubits!(stabilizer::Stabilizer, qubits) # TODO: do we really nead to reset each qubit separately... this is inefficient... can't we just project on all of them at the same time?
#=    s = stabilizer.F22array
    origrows, origcols = size(stabilizer)
    rows = zeros(Bool, origrows)
    for q in qubits
        rows .|= _I.!=s[:,q]
    end
    @assert sum(rows) == length(qubits) "the qubits to be reset are entangled"
    # TODO the check above is not enough
    Stabilizer(stabilizer.phasesF22array[
            .~rows,
            [1,[q+1 for q in 1:origcols if q∉qubits]...]
            ]) # TODO the [qubits...] notation is silly and maybe inefficient
    =#
end


##############################
# Unitary Clifford Operations
##############################

function Base.:(*)(p::AbstractCliffordOperator, s::Stabilizer; phases::Bool=true)
    s = copy(s)
    apply!(s,p; phases=phases)
end

function apply!(s::Stabilizer, p::PauliOperator; phases::Bool=true)
    phases || return s
    for i in 1:length(s.phases)
        s.phases[i] = (s.phases[i]+comm(s[i],p)<<1+p.phase[]<<1)&0x3
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

julia> phase_gate = C\"""Y
                        Z\"""
X ⟼ + Y
Z ⟼ + Z

julia> stab = S\"""XI
                  IZ\""";

julia> entangled = CNOT*stab
+ XX
+ ZZ
```
"""
struct CliffordOperator{Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{UInt64}} <: AbstractCliffordOperator
    phases::Tv
    nqbits::Int
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
    CliffordOperator(vcat((p.phase for p in paulis)...), paulis[1].nqbits, xztox, xztoz)
end

macro C_str(a)
    paulis = [eval(quote @P_str($(strip(s))) end) for s in split(a,'\n')] #TODO seriously!?
    CliffordOperator(paulis)
end

function clifford_transpose(c::CliffordOperator)
    n = c.nqbits
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
    for i in 1:2*c.nqbits
        if i>c.nqbits
            push!(ops,PauliOperator(c.phases[i], ztoxs[:,i-c.nqbits], ztozs[:,i-c.nqbits]))
        else
            push!(ops,PauliOperator(c.phases[i], xtoxs[:,i], xtozs[:,i]))
        end
    end
    ops
end

function Base.getindex(c::CliffordOperator, i::Int)
    xtoxs, ztozs, xtozs, ztoxs = clifford_transpose(c)
    if i>c.nqbits
        PauliOperator(c.phases[i], ztoxs[:,i-c.nqbits], ztozs[:,i-c.nqbits])
    else
        PauliOperator(c.phases[i], xtoxs[:,i], xtozs[:,i])
    end
end

function Base.show(io::IO, c::CliffordOperator)
    xtoxs, ztozs, xtozs, ztoxs = clifford_transpose(c)
    n = c.nqbits
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
    CliffordOperator(copy(c.phases),c.nqbits,copy(c.xztox),copy(c.xztoz))
end

# TODO create Base.permute! and getindex(..., permutation_array)
function permute(c::CliffordOperator,p::AbstractArray{T,1} where T) # TODO this is extremely slow stupid implementation
    ops = getallpaulis_(c)
    CliffordOperator([permute(ops[i],p) for i in 1:2*c.nqbits][vcat(p,p.+c.nqbits)])
end
                  
function apply!(s::Stabilizer, c::CliffordOperator; phases::Bool=true)
    for row_stab in 1:length(s.phases)
        new_stabrowx = zeros(UInt64, length(s.xzs[1,1:end÷2]))
        new_stabrowz = zeros(UInt64, length(s.xzs[1,1:end÷2]))
        for row_clif in 1:s.nqbits
            bigrow = row_clif÷64+1  # TODO use _div and _mod
            smallrow = (row_clif-1)%64  # TODO use _div and _mod
            xztox = c.xztox[row_clif,:] .& s.xzs[row_stab,:]
            xztoz = c.xztoz[row_clif,:] .& s.xzs[row_stab,:]
            new_stabrowx[bigrow] |= xor_bits_(reduce(⊻, xztox)) << smallrow
            new_stabrowz[bigrow] |= xor_bits_(reduce(⊻, xztoz)) << smallrow
            phases && (s.phases[row_stab] = (s.phases[row_stab]+sum(count_zeros, xztoz[1:end÷2] .& xztox[end÷2+1:end])<<1)&0x3)
        end
        s.xzs[row_stab,1:end÷2] .= new_stabrowx
        s.xzs[row_stab,end÷2+1:end] .= new_stabrowz
    end
    s
end

function apply!(s::Stabilizer, c::CliffordOperator, single_qbit_offset::Int)
    bigs = single_qbit_offset÷64+1  # TODO use _div and _mod
    smalls = (single_qbit_offset-1)%64  # TODO use _div and _mod
    lowbit = UInt64(0x1)
    for row_stab in 1:length(s.phases)
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


const CNOT = C"""XX
                 IX
                 ZI
                 ZZ"""

const SWAP = C"""IX
                 XI
                 IZ
                 ZI"""

const Hadamard = C"""Z
                     X"""

const Phase = C"""Y
                  Z"""

const CliffordId = C"""X
                       Z"""

##############################
# Helpers for Clifford Operators
##############################

function (⊗)(l::CliffordOperator, r::CliffordOperator) # TODO this is extremely slow stupid implementation
    opsl = getallpaulis_(l)
    opsr = getallpaulis_(r)
    onel = one(opsl[1])
    oner = one(opsr[1])
    opsl = [l⊗oner for l in opsl]
    opsr = [onel⊗r for r in opsr]
    CliffordOperator(vcat(opsl[1:end÷2],opsr[1:end÷2],opsl[end÷2+1:end],opsr[end÷2+1:end]))
end

function Base.:(*)(l::AbstractCliffordOperator, r::CliffordOperator)
    rstab = Stabilizer(getallpaulis_(r)) # TODO this is a bit awkward... and fragile... turning a CliffordOp into a Stabilizer
    apply!(rstab,l)
    CliffordOperator([rstab[i] for i in 1:length(rstab.phases)])
end

end #module
