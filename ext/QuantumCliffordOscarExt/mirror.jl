struct Mirror{TG, TE} <: AbstractQECC
    G::TG
    A::Vector{TE}
    B::Vector{TE}
    symmetric::Bool

    function Mirror{TG, TE}(
        G::TG, A::Vector{TE}, B::Vector{TE}, symmetric::Bool, ::Val{:validated}
    ) where {TG, TE}
        new{TG, TE}(G, A, B, symmetric)
    end
end

function _mirror_element(G::FinGenAbGroup, coordinates::Tuple{Vararg{Integer}})
    valid = length(coordinates) == length(gens(G)) &&
        all(typemin(Int) <= coordinate <= typemax(Int) for coordinate in coordinates)
    valid || throw(ArgumentError("the coordinates $coordinates do not define an element of G"))
    G(Int[coordinates...])
end

function _mirror_element(G, element)
    element isa typeof(_mirror_identity(G)) && parent(element) === G && return element
    throw(ArgumentError("each member of A and B must be an element of G"))
end

_mirror_identity(G::FinGenAbGroup) = zero(G)
_mirror_identity(G) = one(G)

function _mirror_elements(G, subset)
    elements = typeof(_mirror_identity(G))[]
    for element in subset
        push!(elements, _mirror_element(G, element))
    end
    unique!(elements)
end

function Mirror(G::Union{Group, FinGenAbGroup}, A::AbstractVector, B::AbstractVector; symmetric::Bool=true)
    is_finite(G) || throw(ArgumentError("G must be finite"))
    elements_A = _mirror_elements(G, A)
    elements_B = _mirror_elements(G, B)
    Mirror{typeof(G), typeof(_mirror_identity(G))}(
        G, elements_A, elements_B, symmetric, Val(:validated)
    )
end

Mirror(G, A::AbstractVector, B::AbstractVector, symmetric::Bool) = Mirror(G, A, B; symmetric)

_mirror_product(a::FinGenAbGroupElem, b::FinGenAbGroupElem) = a + b
_mirror_product(a, b) = a * b

_mirror_inverse(a::FinGenAbGroupElem) = -a
_mirror_inverse(a) = inv(a)

function parity_matrix(c::Mirror)
    G, A, B, symmetric = c.G, c.A, c.B, c.symmetric
    elems = collect(G)
    n = length(elems)
    idx = Dict(elems[i] => i for i in 1:n)
    H = falses(n, 2n)

    for (i, g) in enumerate(elems)
        g_inv = _mirror_inverse(g)
        for b in B
            support_x = symmetric ? _mirror_product(b, g_inv) : _mirror_product(g_inv, b)
            H[i, idx[support_x]] = true
        end
        for a in A
            support_z = _mirror_product(a, g)
            H[i, n + idx[support_z]] = true
        end
    end
    H
end
