"""A multiqubit operator corresponding to all identities except for Pauli Z at `i`. See also: [`sY`](@ref), [`sMY`](@ref)"""
function single_z(n,i)
    p = zero(PauliOperator, n)
    p[i] = (false, true)
    p
end

"""A multiqubit operator corresponding to all identities except for Pauli X at `i`. See also: [`sX`](@ref), [`sMX`](@ref)"""
function single_x(n,i)
    p = zero(PauliOperator, n)
    p[i] = (true, false)
    p
end

"""A multiqubit operator corresponding to all identities except for Pauli Y at `i`. See also: [`sY`](@ref), [`sMY`](@ref)"""
function single_y(n,i)
    p = zero(PauliOperator, n)
    p[i] = (true, true)
    p
end

# TODO make faster by using fewer initializations, like in Base.zero
function Base.one(::Type{T}, n; basis=:Z) where {T<:Tableau}# TODO support `basis` in all other `one(::[Mixed][De]Stabilizer)` functions
    if basis==:X
        T(LinearAlgebra.I(n),falses(n,n))
    elseif basis==:Y
        T(LinearAlgebra.I(n),LinearAlgebra.I(n))
    elseif basis==:Z
        T(falses(n,n),LinearAlgebra.I(n))
    else
        throw(ErrorException("`basis` should be one of :X, :Y, or :Z"))
    end
end
Base.one(::Type{<:Stabilizer}, n; basis=:Z) = Stabilizer(one(Tableau,n; basis)) # TODO make it type preserving
Base.one(s::Stabilizer; basis=:Z) = one(Stabilizer, nqubits(s); basis)
Base.one(::Type{<:Destabilizer}, n) = Destabilizer(vcat(one(Tableau, n, basis=:X),one(Tableau, n, basis=:Z)))
function Base.one(::Type{<:MixedStabilizer}, r, n, basis=:Z)
    s = one(Stabilizer, n; basis=basis)
    MixedStabilizer(s,r)
end
function Base.one(::Type{<:MixedDestabilizer}, r, n)
    d = one(Tableau, n; basis=:X)
    s = one(Tableau, n; basis=:Z)
    MixedDestabilizer(vcat(d,s),r)
end
Base.one(T::Type{<:MixedDestabilizer}, n) = one(T, n, n)
function Base.one(c::CliffordOperator)
    n = nqubits(c)
    one(typeof(c),n)
end
Base.one(::Type{<:CliffordOperator}, n) = CliffordOperator(one(Destabilizer,n))

"""Prepare one or more Bell pairs (with optional phases).

```jldoctest
julia> bell()
+ XX
+ ZZ

julia> bell(2)
+ XX__
+ ZZ__
+ __XX
+ __ZZ

julia> bell((true, false))
- XX
+ ZZ

julia> bell([true, false, true, true])
- XX__
+ ZZ__
- __XX
- __ZZ
```
"""
function bell end

function bell()
    copy(S"XX
           ZZ")
end

function bell(n::Int)
    tensor_pow(bell(), n)
end

function bell(phase::Tuple{Bool, Bool})
    s = bell()
    phases(s)[1] = phase[1] ? 0x2 : 0x0
    phases(s)[2] = phase[2] ? 0x2 : 0x0
    s
end

function bell(phases::AbstractVector{Tuple{Bool, Bool}})
    ⊗((bell(t) for t in phases)...)
end

function bell(bellphases::AbstractVector{Bool})
    s = bell(length(bellphases)÷2)
    for i in 1:length(s)
        phases(s)[i] = bellphases[i] ? 0x2 : 0x0
    end
    s
end

"""
Prepare a GHZ state of n qubits.

```jldoctest
julia> ghz()
+ XXX
+ ZZ_
+ _ZZ

julia> ghz(2)
+ XX
+ ZZ

julia> ghz(4)
+ XXXX
+ ZZ__
+ _ZZ_
+ __ZZ
```
"""
function ghz end

function ghz()
    copy(S"XXX
           ZZI
           IZZ")
end

function ghz(n::Int)
    s = zero(Stabilizer, n)
    for i in 1:n
        s[1,i] = (true,false)
    end
    for i in 2:n
        s[i,i] = (false,true)
        s[i,i-1] = (false,true)
    end
    s
end

# TODO document these explicitly
# TODO cluster states, toric code, planar code, other codes from python libraries
