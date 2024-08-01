import Nemo: characteristic, lift, matrix_repr

"""
Classical codes lifted over a permutation group ring [panteleev2021degenerate](@cite) [panteleev2022asymptotically](@cite).

- `A::Matrix`: the base matrix of the code, whose elements are in a permutation group ring.
- `repr::Function`: a function that converts a permutation group ring element to a matrix;
  default to be [`permutation_repr`](@ref) for GF(2)-algebra.

The parity-check matrix is constructed by applying `repr` to each element of `A`, which is mathematically a linear map from a group algebra element to a binary matrix.
This will enlarge the parity check matrix from `A` with each element being inflated into a matrix. The procedure is called a lift [panteleev2022asymptotically](@cite).

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
    !(characteristic(base_ring(A[1,1])) == 2) && error("The default permutation representation applies only to GF(2) group algebra")
    LiftedCode(A, permutation_repr)
end

function LiftedCode(perm_array::Matrix{<:Perm})
    l = length(perm_array[1,1].d)
    R = PermutationGroupRing(GF(2), l)
    A = Matrix(matrix_space(R, size(perm_array)...)(perm_array)) # convert perms to group ring elements
    LiftedCode(A)
end

function LiftedCode(shift_array::Matrix{Int}, l::Int)
    perm_array = map(n->cyclic_permutation(n, l), shift_array)
    LiftedCode(perm_array)
end

"""
Represent a permutation group ring element by mapping permutations to circulant matrices.
"""
function permutation_repr(x::PermGroupRingElem{FqFieldElem})
    mat = zeros(Bool, parent(x).l, parent(x).l)
    for k in keys(x.coeffs)
        c = Int(lift(ZZ, x.coeffs[k])) # the coefficient of the group element `k`, which is converted to Int type for matrix calculation
        mat += c .* Array(matrix_repr(k)) # the coefficient times the matrix representation of the corresponding group element
    end
    return mat
end

function lift(repr::Function, mat::Matrix{PermGroupRingElem})
    vcat([hcat([repr(mat[i, j]) for j in axes(mat, 2)]...) for i in axes(mat, 1)]...)
end

function parity_checks(c::LiftedCode)
    h = lift(c.repr, c.A)
    rk = mod2rank(h) # TODO mod2rank and canonicalize, which is more efficient?
    rk < size(h, 1) && @warn "The lifted code has redundant rows"
    return h
end

code_n(c::LiftedCode) = size(c.A, 2) * size(c.repr(parent(c.A[1,1])(0)), 2)

function mod2rank(h::Matrix{<:Integer}) # TODO move this mod2rank to a common place
    Z2, _ = residue_ring(ZZ, 2)
    S = matrix_space(Z2, size(h)...)
    rank(S(h))
end

code_s(c::LiftedCode) = mod2rank(parity_checks(c))

code_k(c::LiftedCode) = code_n(c) - code_s(c)
