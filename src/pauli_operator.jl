"""
$(TYPEDEF)

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
struct PauliOperator{Tₚ<:AbstractArray{UInt8,0}, Tᵥ<:AbstractVector{<:Unsigned}} <: AbstractCliffordOperator
    phase::Tₚ
    nqubits::Int
    xz::Tᵥ
end

"""$(TYPEDSIGNATURES)"""
PauliOperator(phase::UInt8, nqubits::Int, xz::Tᵥ) where Tᵥ<:AbstractVector{<:Unsigned} = PauliOperator(fill(UInt8(phase),()), nqubits, xz)

"""$(TYPEDSIGNATURES)"""
function PauliOperator(phase::UInt8, x::AbstractVector{Bool}, z::AbstractVector{Bool})
    phase = fill(UInt8(phase),())
    xs = reinterpret(UInt,BitVector(x).chunks)::Vector{UInt}
    zs = reinterpret(UInt,BitVector(z).chunks)::Vector{UInt}
    xzs = cat(xs, zs, dims=1)
    PauliOperator(phase, length(x), xzs)
end

"""$(TYPEDSIGNATURES)"""
PauliOperator(x::AbstractVector{Bool}, z::AbstractVector{Bool}) = PauliOperator(0x0, x, z)

"""$(TYPEDSIGNATURES)"""
PauliOperator(xz::AbstractVector{Bool}) = PauliOperator(0x0, (@view xz[1:end÷2]), (@view xz[end÷2+1:end]))

"""
$(TYPEDSIGNATURES)

Get a view of the X part of the `UInt` array of packed qubits of a given Pauli operator.
"""
function xview(p::PauliOperator)
    @view p.xz[1:end÷2]
end

"""
$(TYPEDSIGNATURES)

Get a view of the Y part of the `UInt` array of packed qubits of a given Pauli operator.
"""
function zview(p::PauliOperator)
    @view p.xz[end÷2+1:end]
end

"""
$(TYPEDSIGNATURES)

Extract as a new bit array the X part of the `UInt` array of packed qubits of a given Pauli operator.
"""
function xbit(p::PauliOperator)
    one = eltype(p.xz)(1)
    size = sizeof(eltype(p.xz))*8
    [(word>>s)&one==one for word in xview(p) for s in 0:size-1][begin:p.nqubits]
end

"""
$(TYPEDSIGNATURES)

Extract as a new bit array the Z part of the `UInt` array of packed qubits of a given Pauli operator.
"""
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

Base.getindex(p::PauliOperator{Tₚ,Tᵥ}, i::Int) where {Tₚ, Tᵥₑ<:Unsigned, Tᵥ<:AbstractVector{Tᵥₑ}} = (p.xz[_div(Tᵥₑ, i-1)+1] & Tᵥₑ(0x1)<<_mod(Tᵥₑ,i-1))!=0x0, (p.xz[end÷2+_div(Tᵥₑ,i-1)+1] & Tᵥₑ(0x1)<<_mod(Tᵥₑ,i-1))!=0x0
Base.getindex(p::PauliOperator{Tₚ,Tᵥ}, r) where {Tₚ, Tᵥₑ<:Unsigned, Tᵥ<:AbstractVector{Tᵥₑ}} = PauliOperator(p.phase[], xbit(p)[r], zbit(p)[r])

function Base.setindex!(p::PauliOperator{Tₚ,Tᵥ}, (x,z)::Tuple{Bool,Bool}, i) where {Tₚ, Tᵥₑ, Tᵥ<:AbstractVector{Tᵥₑ}}
    if x
        p.xz[_div(Tᵥₑ,i-1)+1] |= Tᵥₑ(0x1)<<_mod(Tᵥₑ,i-1)
    else
        p.xz[_div(Tᵥₑ,i-1)+1] &= ~(Tᵥₑ(0x1)<<_mod(Tᵥₑ,i-1))
    end
    if z
        p.xz[end÷2+_div(Tᵥₑ,i-1)+1] |= Tᵥₑ(0x1)<<_mod(Tᵥₑ,i-1)
    else
        p.xz[end÷2+_div(Tᵥₑ,i-1)+1] &= ~(Tᵥₑ(0x1)<<_mod(Tᵥₑ,i-1))
    end
    p
end

Base.firstindex(p::PauliOperator) = 1

Base.lastindex(p::PauliOperator) = p.nqubits

Base.eachindex(p::PauliOperator) = 1:p.nqubits

Base.size(pauli::PauliOperator) = (pauli.nqubits,)

Base.length(pauli::PauliOperator) = pauli.nqubits

nqubits(pauli::PauliOperator) = pauli.nqubits

Base.:(==)(l::PauliOperator, r::PauliOperator) = r.phase==l.phase && r.nqubits==l.nqubits && r.xz==l.xz

Base.hash(p::PauliOperator, h::UInt) = hash(p.phase,hash(p.nqubits,hash(p.xz, h)))

Base.copy(p::PauliOperator) = PauliOperator(copy(p.phase),p.nqubits,copy(p.xz))

"""$(TYPEDSIGNATURES)"""
function LinearAlgebra.inv(p::PauliOperator)
  ph = p.phase[]
  phin = xor((ph << 1) & ~(UInt8(1) << 2), ph)
  return PauliOperator(phin, p.nqubits, copy(p.xz))
end

function Base.deleteat!(p::PauliOperator, subset)
    p =p[setdiff(1:length(p), subset)]
    return p
end

_nchunks(i::Int,T::Type{<:Unsigned}) = 2*( (i-1) ÷ (8*sizeof(T)) + 1 )
Base.zero(::Type{PauliOperator{Tₚ, Tᵥ}}, q) where {Tₚ,T<:Unsigned,Tᵥ<:AbstractVector{T}} = PauliOperator(zeros(UInt8), q, zeros(T, _nchunks(q,T)))
Base.zero(::Type{PauliOperator}, q) = zero(PauliOperator{Array{UInt8, 0}, Vector{UInt}}, q)
Base.zero(p::P) where {P<:PauliOperator} = zero(P, nqubits(p))

"""
$(TYPEDSIGNATURES)

Zero-out the phases and single-qubit operators in a [`PauliOperator`](@ref).
"""
@inline function zero!(p::PauliOperator{Tₚ,Tᵥ}) where {Tₚ, Tᵥₑ<:Unsigned, Tᵥ<:AbstractVector{Tᵥₑ}}
    fill!(p.xz, zero(Tᵥₑ))
    p.phase[] = 0x0
    p
end

"""
$(TYPEDSIGNATURES)

Embed a Pauli operator in a larger Pauli operator.

```jldoctest
julia> embed(5, 3, P"-Y")
- __Y__

julia> embed(5, (3,5), P"-YX")
- __Y_X
```
"""
function embed(n::Int, i::Int, p::PauliOperator)
    if nqubits(p) == 1
        pout = zero(typeof(p), n)
        pout[i] = p[1]
        pout.phase[] = p.phase[]
        return pout
    else
        throw(ArgumentError("""
        You are trying to embed a small Pauli operator into a larger Pauli operator.
        However, you have not given all the positions at which the operator needs to be embedded.
        If you are directly calling `embed`, use the form `embed(nqubits, indices::Tuple, p::PauliOperator)`.
        If you are not using `embed` directly, then `embed` must have been incorrectly called
        by one of the functions you have called.
        """))
    end
end

"""$(TYPEDSIGNATURES)"""
function embed(n::Int, indices, p::PauliOperator)
    if nqubits(p) == length(indices)
        pout = zero(typeof(p), n)
        @inbounds @simd for i in eachindex(indices)
            pout[indices[i]] = p[i]
        end
        pout.phase[] = p.phase[]
        return pout
    else
        throw(ArgumentError(lazy"""
        You are trying to embed a small Pauli operator into a larger Pauli operator.
        However, you have not given all the positions at which the operator needs to be embedded.
        The operator you are embedding is of length $(length(p)), but you have specified $(length(indices)) indices.
        If you are not using `embed` directly, then `embed` must have been incorrectly called
        by one of the functions you have called.
        """))
    end
end
