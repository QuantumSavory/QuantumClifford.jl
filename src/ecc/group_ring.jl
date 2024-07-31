import Base: *, +, -, ==, deepcopy_internal, show, isone, iszero, one, zero, adjoint

import Nemo: Ring, RingElem, RingElement, NCRing, NCRingElem, NCRingElement, Perm, CacheDictType,
    @attributes, base_ring, base_ring_type, elem_type, get_cached!,
    is_domain_type, is_exact_type, parent, parent_type


"""
Permutation [group ring](https://en.wikipedia.org/wiki/Group_ring) over a base ring, also known as group algebra.

- `base_ring`: The base ring, whose elements are used as coefficients of permutations.
- `l`: The length of permutations.

Basic usage:

- To construct a permutation group ring, use [`PermutationGroupRing`](@ref).
- To construct an element of a permutation group ring, use the call syntax of the parent object.
- To perform arithmetic operations, use the standard arithmetic operators.

```jldoctest
julia> ENV["NEMO_PRINT_BANNER"] = "false"; using Nemo: GF, Perm

julia> R = PermutationGroupRing(GF(2), 3)
Permutation group ring over Prime field of characteristic 2

julia> f0 = R(0)
Dict{Perm, Nemo.FqFieldElem}()

julia> f1 = R(1)
Dict{Perm{Int64}, Nemo.FqFieldElem}(() => 1)

julia> f2 = R(Perm([1,3,2]))
Dict{Perm{Int64}, Nemo.FqFieldElem}((2,3) => 1)

julia> f3 = R(Dict(Perm(3) => 1, Perm([3,1,2]) => 1))
Dict{Perm{Int64}, Nemo.FqFieldElem}(() => 1, (1,3,2) => 1)

julia> f1 + f2
Dict{Perm, Nemo.FqFieldElem}(() => 1, (2,3) => 1)

julia> f2 * f3
Dict{Perm, Nemo.FqFieldElem}((1,3) => 1, (2,3) => 1)

julia> f3 * 1 + 1
Dict{Perm, Nemo.FqFieldElem}((1,3,2) => 1)

julia> f3'
Dict{Perm, Nemo.FqFieldElem}(() => 1, (1,2,3) => 1)
```

See also: [`PermGroupRingElem`](@ref).
"""
@attributes mutable struct PermGroupRing{T<:RingElement} <: NCRing
    base_ring::Ring
    l::Int

    function PermGroupRing{T}(R::Ring, l::Int, cached::Bool) where {T<:RingElement}
        return get_cached!(PermGroupRingElemID, (R, l), cached) do
            new{T}(R, l)
        end::PermGroupRing{T}
    end
end

const PermGroupRingElemID = CacheDictType{Tuple{Ring, Int}, NCRing}()

"""
Element of a [`PermGroupRing`](@ref).

- `coeffs`: A dictionary of permutations and their coefficients. Empty dictionary represents zero.
- `parent`: The parent group ring in type `PermGroupRing`.

See also: [`PermGroupRing`](@ref).
"""
mutable struct PermGroupRingElem{T<:RingElement} <: NCRingElem
    coeffs::Dict{<:Perm,T}
    parent::PermGroupRing{T}

    function PermGroupRingElem{T}(coeffs::Dict{<:Perm,T}) where {T<:RingElement}
        filter!(x -> !iszero(x[2]) , coeffs) # remove zeros
        return new{T}(coeffs)
    end

    function PermGroupRingElem{T}() where {T<:RingElement}
        return new{T}(Dict{Perm,T}())
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

# String I/O

function show(io::IO, R::PermGroupRing)
    print(io, "Permutation group ring over ")
    show(io, base_ring(R))
 end

 function show(io::IO, f::PermGroupRingElem)
    print(io, f.coeffs)
 end

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

# TODO Some functionality are expected by ring interfaces but not necessary for ECC construction,
# including hash, exact division, random generation, promotion rules

# Constructors by overloading the call syntax for parent objects

function (R::PermGroupRing{T})() where {T<:RingElement}
    r = PermGroupRingElem{T}()
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(coeffs::Dict{<:Perm,<:Union{T,Integer,Rational,AbstractFloat}}) where {T<:RingElement}
    if valtype(coeffs) == T
        for (k,v) in coeffs
            length(k.d) == R.l || error("Invalid permutation length")
            parent(v) == R || error("Unable to coerce a group ring element")
        end
    else
        coeffs = Dict(k => base_ring(R)(v) for (k, v) in coeffs)
    end
    r = PermGroupRingElem{T}(coeffs)
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(n::Union{Integer,Rational,AbstractFloat}) where {T<:RingElement}
    coeffs = iszero(n) ? Dict{Perm,T}() : Dict(Perm(R.l) => base_ring(R)(n))
    r = PermGroupRingElem{T}(coeffs)
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(n::T) where {T<:RingElement}
    base_ring(R) == parent(n) || error("Unable to coerce a group ring element")
    r = PermGroupRingElem{T}(Dict(Perm(R.l) => n))
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(p::Perm) where {T<:RingElement}
    length(p.d) == R.l || error("Invalid permutation length")
    r = PermGroupRingElem{T}(Dict(p => one(base_ring(R))))
    r.parent = R
    return r
end

function (R::PermGroupRing{T})(f::PermGroupRingElem{T}) where {T<:RingElement}
    parent(f) == R || error("Unable to coerce a group ring element")
    return f
end

"""
Permutation group ring constructor.

See also: [`PermutationGroupRing`](@ref).
"""
function PermutationGroupRing(R::Ring, l::Int, cached::Bool=true)
    T = elem_type(R)
    return PermGroupRing{T}(R, l, cached)
end

# adjoint

function adjoint(a::PermGroupRingElem{T}) where {T<:FqFieldElem}
    r = parent(a)()
    for (k, v) in a.coeffs
        r.coeffs[inv(k)] = v
    end
    return r
end

"""
Construct a cyclic permutation of length `l` with a shift `n`.
"""
cyclic_permutation(n::Int, l::Int) = Perm(vcat(n+1:l, 1:n))
