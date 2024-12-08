"""Hypergraph product of two classical codes."""
function hgp(h₁,h₂)
    r₁, n₁ = size(h₁)
    r₂, n₂ = size(h₂)
    hx = hcat(kron(h₁, LinearAlgebra.I(n₂)), kron(LinearAlgebra.I(r₁), h₂'))
    hz = hcat(kron(LinearAlgebra.I(n₁), h₂), kron(h₁', LinearAlgebra.I(r₂)))
    hx, hz
end

"""Implemented in a package extension with Hecke."""
function check_repr_commutation_relation end

"""Implemented in a package extension with Hecke."""
function check_commutative_group_algebra end
