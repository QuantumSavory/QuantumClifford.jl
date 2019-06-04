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

export @P_str, PauliOperator, ⊗, I, X, Y, Z,
    @S_str, Stabilizer, prodphase, comm, ⊕, isfullstabilizer, canonicalize!,
    generate!, project!, reset_qubits!, traceout_qubits!,
    apply!,
    CliffordOperator, @C_str, CNOT, SWAP, Hadamard, Phase, CliffordId

# Predefined constants representing the single qubit Pauli operators encoded
# as elements of ``F(2,2)`` in the low bits of UInt8.
const _I = 0b00
const _X = 0b10
const _Y = 0b11
const _Z = 0b01

# Predefined constants representing the permitted phases encoded
# as elements of ``F(2)`` in the low bits of UInt8.
const _p  = 0x00
const _pi = 0x01
const _m  = 0x02
const _mi = 0x03

# Dictionaries used for parsing and printing.
const l2F = Dict('I'=>_I,'X'=>_X,'Z'=>_Z,'Y'=>_Y)
const F2l = Dict(v=>k for (k,v) in l2F)
F2l[_I] = '_'
l2F['_'] = _I
const phasedict = Dict(""=>_p,"+"=>_p,"i"=>_pi,"+i"=>_pi,"-"=>_m,"-i"=>_mi)

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

Internally, each single-qubit operator is an element of `F(2,2)` encoded
in the low bits of a `UInt8`. The phase is encoded as element of `F(4)`, also
in a `UInt8`. The implementation uses a single array with properties that
provide views into the array. The current implementation wastes quite a bit of
space.

```jldoctest
julia> p = P"-IZXY";

julia> p.phase, p.F22array, p.phaseF22array
(0x02, UInt8[0x00, 0x01, 0x02, 0x03], UInt8[0x02, 0x00, 0x01, 0x02, 0x03])
```
"""
struct PauliOperator{T<:AbstractArray{UInt8,1}} <: AbstractCliffordOperator
    phaseF22array::T # the first element is the phase (0,1,2,3 for +,+i,-,-i)
end

macro P_str(a)
    f22array = collect(l2F[l] for l in filter(x->occursin(x,"IZXY"),a))
    phase = phasedict[strip(filter(x->!occursin(x,"IZXY"),a))]
    PauliOperator(vcat([phase],f22array))
end

function Base.getproperty(pauli::PauliOperator, name::Symbol)
    if name==:phase
        pauli.phaseF22array[1]
    elseif name==:F22array
        @view pauli.phaseF22array[2:end]
    else
        getfield(pauli, name)
    end
end

function Base.setproperty!(pauli::PauliOperator, name::Symbol, x)
    if name==:phase
        pauli.phaseF22array[1] = x
    elseif name==:F22array
        pauli.phaseF22array[2:end] = x
    else
        setfield!(pauli, name, x)
    end
end

Base.propertynames(pauli::PauliOperator, private=false) = (:phase,:F22array,:phaseF22array)

Base.size(pauli::PauliOperator) = size(pauli.F22array)[1]

array2str(a) = join(F2l[l] for l in a)

Base.show(io::IO, p::PauliOperator) = print(io, ["+ ","+i","- ","-i"][p.phase+1]*array2str(p.F22array))

Base.:(==)(l::PauliOperator, r::PauliOperator) = r.phaseF22array==l.phaseF22array

Base.hash(p::PauliOperator, h::UInt) = hash(p.phaseF22array, h)

Base.copy(p::PauliOperator) = PauliOperator(copy(p.phaseF22array))

Base.one(p::PauliOperator) = PauliOperator(zeros(UInt8,size(p)+1))

##############################
# Pauli Operator Helpers
##############################

function prodphase_(l::UInt8, r::UInt8)::UInt8 # TODO this really needs a neater non-lookup representation
    multiplication_table_s = Array{UInt8,2}(   # related TODO - frequently we care only about +/- phases, where a non-lookup is easier
    # I Z X Y  <-- second argument
    [[0 0 0 0]; # I  <-- first argument
     [0 0 1 3]; # Z
     [0 3 0 1]; # X
     [0 1 3 0]] # Y
    )
    multiplication_table_s[l+1,r+1]
end

prodphase_(l::AbstractArray{UInt8,1}, r::AbstractArray{UInt8,1})::UInt8 = sum(prodphase_.(l,r))

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
prodphase(l::PauliOperator, r::PauliOperator)::UInt8 = (l.phase+r.phase+prodphase_(l.F22array,r.F22array))%4

function comm_(l::UInt8, r::UInt8)::UInt8 # based on the twisted product from arxiv 0304161
    (l&_Z)&((r&_X)>>1) + (r&_Z)&((l&_X)>>1)
end

comm_(l::PauliOperator, r::PauliOperator)::UInt8 = sum(comm_.(l.F22array,r.F22array))

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
comm(l, r)::UInt8 = comm_(l,r)%2


function Base.:(*)(l::PauliOperator, r::PauliOperator)
    p = copy(l)
    p.F22array .⊻= r.F22array
    p.phase = prodphase(l,r)
    p
end

(⊗)(l::PauliOperator, r::PauliOperator) = PauliOperator(vcat([(l.phase+r.phase)%0x4], l.F22array, r.F22array))

function Base.:(*)(l, r::PauliOperator)
    p = copy(r)
    if l==1
        nothing
    elseif l==1im
        p.phase = (p.phase + 1)%4
    elseif l==-1
        p.phase = (p.phase + 1)%4
    elseif l==-1im
        p.phase = (p.phase + 1)%4
    else
        throw(DomainError(l,"Only {±1,±i} are permitted as phases."))
    end
    p
end

Base.:(+)(p::PauliOperator) = p

function Base.:(-)(p::PauliOperator)
    p = copy(p)
    p.phase = (p.phase+2)%4
    p
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
It encodes the stabilizer in a 2D array of single-qubit operators encoded
as elements of ``F(2,2)``. Rows correspond to multi-qubit Pauli operators and
columns correspond to qubits. Instances can be created with the `S""` literal.

There are no automatic checks for correctness (i.e. independence of all rows,
commutativity of all rows, hermiticity of all rows).

See also: [`PauliOperator`](@ref), [`canonicalize!`](@ref)
"""
struct Stabilizer{T<:AbstractArray{UInt8,2}}
    phasesF22array::T
end

function parse_listofpaulis(a)::Array{UInt8,2}
    f22array = hcat((collect(l2F[l] for l in filter(x->occursin(x,"_IZXY"),s)) for s in split(a,'\n'))...)'
    phases = collect(phasedict[strip(filter(x->!occursin(x,"_IZXY"),l))] for l in split(a,'\n'))
    hcat(phases,f22array)
end

macro S_str(a)
    Stabilizer(parse_listofpaulis(a))
end

function Base.getproperty(stab::Stabilizer, name::Symbol)
    if name==:phases
        @view stab.phasesF22array[:,1]
    elseif name==:F22array
        @view stab.phasesF22array[:,2:end]
    else
        getfield(stab, name)
    end
end

function Base.setproperty!(stab::Stabilizer, name::Symbol, x)
    if name==:phases
        stab.phasesF22array[:,1] = x
    elseif name==:F22array
        stab.phasesF22array[:,2:end] = x
    else
        setfield!(stab, name, x)
    end
end

Base.propertynames(stab::Stabilizer, private=false) = (:phases,:F22array,:phasesF22array)

Base.getindex(stab::Stabilizer, i::Int) = PauliOperator(@view stab.phasesF22array[i,:])
Base.getindex(stab::Stabilizer, r::UnitRange) = Stabilizer(@view stab.phasesF22array[r,:])

function Base.setindex!(stab::Stabilizer, pauli::PauliOperator, i)
    stab.phasesF22array[i,:] = pauli.phaseF22array
    pauli
end

Base.firstindex(stab::Stabilizer) = 1

Base.lastindex(stab::Stabilizer) = size(stab.phasesF22array,1)

Base.size(stabilizer::Stabilizer, args...) = size(stabilizer.F22array, args...)

Base.show(io::IO, s::Stabilizer) = print(io,
                                         join(map(t->t[1]*t[2],
                                                  zip(map(b->["+ ","+i","- ","-i"][b+1], s.phases),
                                                      mapslices(array2str, s.F22array, dims=2))),
                                              '\n'))

Base.:(==)(l::Stabilizer, r::Stabilizer) = r.phasesF22array==l.phasesF22array

Base.hash(s::Stabilizer, h::UInt) = hash(s.phasesF22array, h)

Base.copy(s::Stabilizer) = Stabilizer(copy(s.phasesF22array))

function rowswap!(s::Stabilizer, i, j) # Written only so we can avoid copying in `canonicalize!`
    if i == j
        return
    end
    for k in 1:size(s.phasesF22array,2)
        s.phasesF22array[i,k], s.phasesF22array[j,k] = s.phasesF22array[j,k], s.phasesF22array[i,k]
    end
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
function canonicalize!(stabilizer::Stabilizer) # TODO simplify by using the new Pauli like interface instead of f22array
    f22array = stabilizer.F22array  
    rows, columns = size(stabilizer)
    i = 1
    for j in 1:columns
        k = findfirst(e->e&_X!=0, f22array[i:end,j]) # if X or Y
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i)
            for m in 1:rows
                if f22array[m,j]&_X!=0 && m!=i # if X or Y
                    stabilizer[m] = stabilizer[m] * stabilizer[i] # TODO this should be in-place
                end
            end
            i += 1
        end
    end
    for j in 1:columns
        k = findfirst(e->e==_Z, f22array[i:end,j]) # if Z
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i)
            for m in 1:rows
                if f22array[m,j]&_Z!=0 && m!=i # if Z or Y
                    stabilizer[m] = stabilizer[m] * stabilizer[i] # TODO this should be in-place
                end
            end
            i += 1
        end
    end
    stabilizer
end

function ishermitian() # TODO write it both for paulis and stabilizers
end

function isfullstabilizer(stabilizer::Stabilizer)
    s = stabilizer.F22array
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
    return true
end

function ⊕(l::Stabilizer, r::Stabilizer)
    stabs = cat(l.F22array,r.F22array,dims=(1,2))
    phases = vcat(l.phases,r.phases)
    Stabilizer(hcat(phases,stabs))
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
function generate!(pauli::PauliOperator, stabilizer::Stabilizer)
    rows, columns = size(stabilizer)
    p = pauli.F22array
    s = stabilizer.F22array
    used_indices = Int[]
    used = 0
    # remove Xs
    while (i=findfirst(e->(e&_X)!=0,p)) !== nothing
        used += findfirst(e->(e&_X)!=0,s[used+1:end,i])
        pauli.phaseF22array .= (pauli*stabilizer[used]).phaseF22array # TODO, this is just a silly way to write it... learn more about broadcast
        push!(used_indices, used)
    end
    while (i=findfirst(e->e==_Z,p)) !== nothing
        used += findfirst(e->e==_Z,s[used+1:end,i])
        pauli.phaseF22array .= (pauli*stabilizer[used]).phaseF22array # TODO, this is just a silly way to write it... learn more about broadcast
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
function project!(stabilizer::Stabilizer,pauli::PauliOperator;keep_result=true)
    # @assert !all(pauli.F22array .== I)
    # @assert pauli.phase ∈ [0x0, 0x2]
    anticommutes = 0                                           
    n = size(stabilizer)[1]
    for i in 1:n
        if comm(pauli,stabilizer[i])!=0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        if keep_result
            canonicalize!(stabilizer)
            new_pauli, _ = generate!(copy(pauli), stabilizer)
            result = new_pauli.phase
        else
            result = nothing
        end
    else
        for i in anticommutes+1:n
            if comm(pauli,stabilizer[i])!=0
                stabilizer[i] = stabilizer[i] * stabilizer[anticommutes] # TODO this should be in-place
                # @assert comm(pauli,stabilizer[i])==0 # this assert is always true in the absence of bugs
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
"""=# #TODO
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

"""
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
""" # TODO it is not exactly mutable, rather it mutates and then returns a view...
function traceout_qubits!(stabilizer::Stabilizer, qubits) # TODO: do we really nead to reset each qubit separately... this is inefficient... can't we just project on all of them at the same time?
    s = stabilizer.F22array
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
end


##############################
# Unitary Clifford Operations
##############################

function Base.:(*)(p::AbstractCliffordOperator, s::Stabilizer)
    s = copy(s)
    apply!(s,p)
end

function apply!(s::Stabilizer, p::PauliOperator)
    for i in 1:size(s)[1]
        s.phases[i] = (s.phases[i] + 2*p.phase + 2*sum((s.F22array[i,:] .!= p.F22array) .& (_I .!= p.F22array) .& (s.F22array[i,:] .!= _I)))%4
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

julia> apply!(stab,phase_gate,(1,)) # if the gate is smaller than the stabilizer, specify the qubits to act on
+ Y_
+ _Z
```
"""
struct CliffordOperator{T<:AbstractArray{UInt8,2}} <: AbstractCliffordOperator
    phasesF22array::T
end

macro C_str(a)
    CliffordOperator(parse_listofpaulis(a))
end

function Base.show(io::IO, c::CliffordOperator)
    a = c.phasesF22array
    n = size(a,2)-1
    for i in 1:n
        print(io, repeat("_",i-1),"X",repeat("_",n-i), " ⟼ ")
        print(io, ["+ ","+i","- ","-i"][a[i,1]+1])
        print(io, array2str(a[i,2:end]))
        println(io)
    end
    for i in 1:n
        print(io, repeat("_",i-1),"Z",repeat("_",n-i), " ⟼ ")
        print(io, ["+ ","+i","- ","-i"][a[i+n,1]+1])
        print(io, array2str(a[i+n,2:end]))
        println(io)
    end
end

function Base.copy(c::CliffordOperator)
    CliffordOperator(copy(c.phasesF22array))
end

function Base.permute!(c::CliffordOperator,p::AbstractArray{T,1} where T)
    nbops = size(c.phasesF22array,1)÷2
    c.phasesF22array .= c.phasesF22array[[p...,(p .+ nbops)...],[1,(1 .+ p)...]]
    c
end

Base.getindex(stab::CliffordOperator, i::Int) = PauliOperator(@view stab.phasesF22array[i,:])                        
                        
function apply!(s::Stabilizer, c::CliffordOperator) # TODO Is this writen in an incredibly inefficient manner? Why is it not just a matrix multiplicaiton?
    rows,cols = size(s)
    for row in 1:rows
        tmp = one(s[1])
        orig_phase = s.phases[row]
        for col in 1:cols
            if s.F22array[row,col] == _X
                tmp.phaseF22array .= (tmp*c[col]).phaseF22array
                # TODO, this above is just a silly way to write it... learn broadcasting
            elseif s.F22array[row,col] == _Z
                tmp.phaseF22array .= (tmp*c[col+cols]).phaseF22array
                # TODO, this above is just a silly way to write it... learn broadcasting
            elseif s.F22array[row,col] == _Y
                tmp.phaseF22array .= (1im*tmp*c[col]*c[col+cols]).phaseF22array
                # TODO, this above is just a silly way to write it... learn broadcasting
            end
        end
        s[row] = tmp
        s.phases[row] = (s.phases[row]+orig_phase)%0x4 # TODO is this the cleanest API
    end
    s
end

function apply!(s::Stabilizer, c::AbstractCliffordOperator, qubits) # TODO qubits should be typed into whatever abstract indexing type exits
    indices = [1, (qubits .+ 1)...]
    view_s = Stabilizer(@view s.phasesF22array[:,indices])
    # TODO overload s[row_indexer,col_indexer] and @view s[row_indexer,col_indexer] so that the above line is simpler
    apply!(view_s, c)
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

function (⊗)(l::CliffordOperator, r::CliffordOperator)
    l,r = l.phasesF22array,r.phasesF22array
    phases = vcat(l[1:end÷2,1], r[1:end÷2,1], l[end÷2+1:end,1], r[end÷2+1:end,1])
    lx,lz = l[1:end÷2,2:end], l[end÷2+1:end,2:end]
    rx,rz = r[1:end÷2,2:end], r[end÷2+1:end,2:end]
    ops = vcat(cat(lx,rx,dims=(1,2)),cat(lz,rz,dims=(1,2)))
    CliffordOperator(hcat(phases,ops))
end

function Base.:(*)(l::AbstractCliffordOperator, r::CliffordOperator)
    r = Stabilizer(copy(r.phasesF22array)) # TODO this is a bit awkward... and fragile... turning a CliffordOp into a Stabilizer
    apply!(r,l)
    CliffordOperator(r.phasesF22array)
end

end #module





















