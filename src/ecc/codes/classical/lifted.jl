import Nemo: characteristic, lift, matrix_repr

"""
Classical codes lifted over a permutation group ring [panteleev2021degenerate](@cite) [panteleev2022asymptotically](@cite).

- `A::Matrix`: the base matrix of the code, whose elements are in a permutation group ring.
- `repr::Function`: a function that converts a permutation group ring element to a matrix;
  default to be [`permutation_repr`](@ref) for GF(2)-algebra.

The parity-check matrix is constructed by applying `repr` to elements of `A`, known as a lift.

See also: [`LPCode`](@ref), [`PermGroupRing`](@ref).
"""
struct LiftedCode <: ClassicalCode
    A::Matrix{PermGroupRingElem}
    repr::Function

    function LiftedCode(A::Matrix{PermGroupRingElem{T}}, repr::Function) where T
        new(A, repr)
    end
end

function LiftedCode(A::Matrix{PermGroupRingElem{FqFieldElem}})
    # TODO we may also want to handle the case where the base ring is not GF(2),
    # such as residue_ring(ZZ, 2), which is mathematically equivalent, and residue_ring(ZZ, n) for n > 2
    !(characteristic(base_ring(A[1,1])) == 2) && error("The default permutation representation applies only to GF(2) group algebra")
    LiftedCode(A, permutation_repr)
end

"""
Represent a permutation group ring element by mapping permutations to circulant matrices.
"""
function permutation_repr(x::PermGroupRingElem{FqFieldElem})
    return sum([Int(lift(ZZ, x.coeffs[k])) .* Array(matrix_repr(k)) for k in keys(x.coeffs)], init=zeros(Bool, parent(x).l, parent(x).l))
end

function lift(repr::Function, mat::Matrix{PermGroupRingElem})
    vcat([hcat([repr(mat[i, j]) for j in axes(mat, 2)]...) for i in axes(mat, 1)]...)
end

parity_checks(c::LiftedCode) = lift(c.repr, c.A)

code_n(c::LiftedCode) = size(c.A, 2) * size(c.repr(parent(c.A[1,1])(0)), 2)

function mod2rank(h::Matrix{<:Integer})
    Z2, _ = residue_ring(ZZ, 2)
    S = matrix_space(Z2, size(h)...)
    rank(S(h))
end

code_s(c::LiftedCode) = mod2rank(parity_checks(c)) # note that redundant rows exist in general

code_k(c::LiftedCode) = code_n(c) - code_s(c)
