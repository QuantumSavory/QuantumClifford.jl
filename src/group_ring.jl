import AbstractAlgebra: Ring, RingElem, add!, addeq!, base_ring, base_ring_type,
    canonical_unit, characteristic, divexact, elem_type, expressify, get_cached!,
    is_domain_type, is_exact_type, is_unit, isequal, mul!, parent, parent_type, zero!
import Base:
    *, +, -, ==, ^, deepcopy_internal, hash, inv, isone, iszero, one, rand, show, zero

using AbstractAlgebra
using Random: Random, GLOBAL_RNG, SamplerTrivial
using RandomExtensions: RandomExtensions, AbstractRNG, Make2

@attributes mutable struct PermGroupRing{T<:RingElement} <: NCRing
    base_ring::Ring
    l::Int

    function PermGroupRing{T}(R::Ring, l::Int, cached::Bool) where {T<:RingElement}
        return get_cached!(PermGroupRingElemID, R, cached) do
            new{T}(R, l)
        end::PermGroupRing{T}
    end
end

const PermGroupRingElemID = AbstractAlgebra.CacheDictType{NCRing,PermGroupRing}()

mutable struct PermGroupRingElem{T<:RingElement} <: NCRingElem
    coeffs::Dict{<:Perm,T}
    parent::PermGroupRing{T}

    function PermGroupRingElem{T}(coeffs::Dict{<:Perm,T}) where {T<:RingElement}
        filter!(x -> x[2] != 0, coeffs)
        return new{T}(coeffs)
    end

    function PermGroupRingElem{T}() where {T<:RingElement}
        return new{T}(Dict{Perm,T}())
    end

    function PermGroupRingElem{T}(n::T, l::Int) where {T<:RingElement}
        if iszero(n)
            coeffs = Dict{Perm,T}()
        else
            coeffs = Dict(Perm(l) => n)
        end
        return new{T}(coeffs)
    end

    function PermGroupRingElem{T}(p::Perm, base_ring::Ring) where {T<:RingElement}
        return new{T}(Dict(p => one(base_ring)))
    end
end

# Data type and parent object methods

parent_type(::Type{PermGroupRingElem{T}}) where {T<:RingElement} = PermGroupRing{T}

elem_type(::Type{PermGroupRing{T}}) where {T<:RingElement} = PermGroupRingElem{T}

base_ring_type(::Type{PermGroupRing{T}}) where {T<:RingElement} = parent_type(T)

base_ring(R::PermGroupRing) = R.base_ring::base_ring_type(R)

parent(f::PermGroupRingElem) = f.parent

is_domain_type(::Type{PermGroupRingElem{T}}) where {T<:RingElement} = is_domain_type(T)

is_exact_type(::Type{PermGroupRingElem{T}}) where {T<:RingElement} = is_exact_type(T)

function hash(f::PermGroupRingElem, h::UInt)
    r = 0x65125ab8e0cd44ca # TODO: what to do with this?
    return xor(r, hash(f.c, h))
end

function deepcopy_internal(f::PermGroupRingElem{T}, dict::IdDict) where {T<:RingElement}
    r = PermGroupRingElem{T}(deepcopy_internal(f.coeffs, dict))
    r.parent = f.parent # parent should not be deepcopied
    return r
end

# Basic manipulation

zero(R::PermGroupRing) = R()

one(R::PermGroupRing) = R(1)

iszero(f::PermGroupRingElem) = f == 0

isone(f::PermGroupRingElem) = f == 1

# Arithmetic functions

function -(a::PermGroupRingElem{T}) where {T<:RingElement}
    r = parent(a)()
    for (k, v) in a.coeffs
        r.coeffs[k] = -v
    end
    return r
end

function +(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement}
    r = parent(a)()
    for (k, v) in a.coeffs
        r.coeffs[k] = v
    end
    for (k, v) in b.coeffs
        if haskey(r.coeffs, k)
            r.coeffs[k] += v
            if r.coeffs[k] == 0
                delete!(r.coeffs, k)
            end
        else
            r.coeffs[k] = v
        end
    end
    return r
end

-(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement} = a + (-b)

function *(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement}
    r = parent(a)()
    for (k1, v1) in a.coeffs
        for (k2, v2) in b.coeffs
            k = k1 * k2
            if haskey(r.coeffs, k)
                r.coeffs[k] += v1 * v2
            else
                r.coeffs[k] = v1 * v2
            end
        end
    end
    filter!(x -> x[2] != 0, r.coeffs)
    return r
end

# Ad hoc arithmetic functions

*(a::PermGroupRingElem{T}, n::T) where {T<:RingElement} = a * parent(a)(n)

*(n::T, a::PermGroupRingElem{T}) where {T<:RingElement} = a * n

+(a::PermGroupRingElem{T}, n::T) where {T<:RingElement} = a + parent(a)(n)

+(n::T, a::PermGroupRingElem{T}) where {T<:RingElement} = a + n

*(a::PermGroupRingElem{T}, p::Perm) where {T<:RingElement} = a * parent(a)(p)

*(p::Perm, a::PermGroupRingElem{T}) where {T<:RingElement} = parent(a)(p) * a

+(a::PermGroupRingElem{T}, p::Perm) where {T<:RingElement} = a + parent(a)(p)

+(p::Perm, a::PermGroupRingElem{T}) where {T<:RingElement} = a + p

# Comparison

function ==(a::PermGroupRingElem{T}, b::PermGroupRingElem{T}) where {T<:RingElement}
    if length(a.coeffs) != length(b.coeffs)
        return false
    end
    for (k, v) in a.coeffs
        if !haskey(b.coeffs, k) || b.coeffs[k] != v
            return false
        end
    end
    return true
end

# Ad hoc comparison

==(a::PermGroupRingElem{T}, n::T) where {T<:RingElement} = length(a.coeffs) == 1 && a.coeffs[Perm(parent(a).l)] == n

==(n::T, a::PermGroupRingElem{T}) where {T<:RingElement} = a == n

==(a::PermGroupRingElem{T}, p::Perm) where {T<:RingElement} = length(a.coeffs) == 1 && a.coeffs[p] == base_ring(parent(a))(1)

==(p::Perm, a::PermGroupRingElem{T}) where {T<:RingElement} = a == p


# TODO Some ring interfaces, which are not required in code construction
# exact division, random generation

# TODO Promotion rules

# Constructors by overloading the call syntax for parent objects

function (R::PermGroupRing{T})() where {T<:RingElement}
    r = PermGroupRingElem{T}()
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(coeffs::Dict{<:Perm,T}) where {T<:RingElement}
    r = PermGroupRingElem{T}(coeffs)
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(coeffs::Dict{<:Perm,<:Union{Integer,Rational,AbstractFloat}}) where {T<:RingElement}
    r = PermGroupRingElem{T}(Dict(k => base_ring(R)(v) for (k, v) in coeffs))
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(n::Union{Integer,Rational,AbstractFloat}) where {T<:RingElement}
    r = PermGroupRingElem{T}(base_ring(R)(n), R.l)
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(n::T) where {T<:RingElement}
    base_ring(R) != parent(n) && error("Unable to coerce group ring element")
    r = PermGroupRingElem{T}(n, l)
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(p::Perm) where {T<:RingElement}
    r = PermGroupRingElem{T}(p, base_ring(R))
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(f::PermGroupRingElem{T}) where {T<:RingElement}
    parent(f) != R || error("Unable to coerce group ring")
    return f
end

# TODO We may need more constructors to remove ambiguities

# Parent constructor

function PermutationGroupRing(R::Ring, l::Int, cached::Bool=true)
    T = elem_type(R)
    return PermGroupRing{T}(R, l, cached)
end
