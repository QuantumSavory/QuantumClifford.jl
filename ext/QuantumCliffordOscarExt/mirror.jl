"""
    Mirror(G, A, B; symmetric=true)

$TYPEDEF

A mirror quantum LDPC code constructed from a finite group `G` and two subsets
`A` and `B`, as introduced in [khesin2026mirror](@cite).

The qubits and stabilizer checks are both indexed by the elements `g` of `G`.
For a symmetric mirror code, the binary support of the check at `g` is

```math
S(g) = Z(A g) X(B g^{-1}).
```

For an asymmetric mirror code (`symmetric=false`), it is

```math
S(g) = Z(A g) X(g^{-1} B).
```

Thus, each check has weight at most `length(A) + length(B)`. The two
constructions coincide for Abelian groups. For a non-Abelian group, the
constructor verifies that all checks commute and throws an `ArgumentError` if
`A` and `B` do not define a stabilizer code. [`parity_checks`](@ref) selects
consistent signs for redundant generators, so their group does not contain
negative identity.

Elements of `A` and `B` can be group elements. For a `FinGenAbGroup`, integer
coordinate tuples are also accepted. Duplicate elements are discarded, and
either subset can be empty.

[`parity_matrix`](@ref) returns one row per group element and `2 * length(G)`
Boolean columns in `(X | Z)` order.

Here is the `[[16, 4, 4]]` code from Table I of
[khesin2026mirror](@cite):

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> G = abelian_group([16]);

julia> c = Mirror(G, [(0,), (4,)], [(1,), (3,), (5,), (11,)]);

julia> code_n(c), code_k(c)
(16, 4)
```

### Fields

$TYPEDFIELDS
"""
struct Mirror{TG, TE} <: AbstractQECC
    """The finite group that indexes the qubits and checks."""
    G::TG
    """The subset that specifies the `Z` support."""
    A::Vector{TE}
    """The subset that specifies the `X` support."""
    B::Vector{TE}
    """Whether to use the symmetric rather than asymmetric construction."""
    symmetric::Bool
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
    code = Mirror{typeof(G), typeof(_mirror_identity(G))}(G, elements_A, elements_B, symmetric)

    if !is_abelian(G) && !QuantumClifford.check_allrowscommute(Stabilizer(parity_matrix(code)))
        throw(ArgumentError("A and B do not define commuting mirror checks for G"))
    end

    code
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

function parity_checks(c::Mirror)
    checks = Stabilizer(parity_matrix(c))
    # Independent canonical rows cannot generate -I. Match the sign of every
    # original support to the product of canonical rows that generates it.
    canonical, _, rank = canonicalize!(copy(checks), ranks=true)
    basis = canonical[1:rank]

    for i in eachindex(checks)
        remainder = generate!(copy(checks[i]), basis; saveindices=false)
        @assert !isnothing(remainder)
        phases(checks)[i] = remainder.phase[]
    end

    checks
end
