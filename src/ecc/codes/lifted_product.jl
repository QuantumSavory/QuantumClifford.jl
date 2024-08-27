import Nemo: matrix_space
import LinearAlgebra

"""
Lifted product codes [panteleev2021degenerate](@cite) [panteleev2022asymptotically](@cite)

- `A::Matrix{PermGroupRingElem}`: the first base matrix for constructing the lifted product code, whose elements are in a permutation group ring;
- `B::Matrix{PermGroupRingElem}`: the second base matrix for constructing the lifted product code, whose elements are in the same permutation group ring as `A`;
- `repr::Function`: a function that converts the permutation group ring element to a matrix;
 default to be [`permutation_repr`](@ref) for GF(2)-algebra.

A lifted product code is constructed by hypergraph product of the two lifted codes `c₁` and `c₂`.
Here, the hypergraph product is taken over a group ring, which serves as the base ring for both lifted codes.
After the hypergraph product, the parity-check matrices are lifted by `repr`.
The lifting is achieved by applying `repr` to each element of the matrix resulted from the hypergraph product, which is mathematically a linear map from a group algebra element to a binary matrix.

See also: [`LiftedCode`](@ref), [`PermGroupRing`](@ref).
"""

struct LPCode <: AbstractECC
    A::PermGroupRingMatrix
    B::PermGroupRingMatrix
    repr::Function

    function LPCode(A::PermGroupRingMatrix, B::PermGroupRingMatrix, repr::Function)
        A[1, 1].parent == B[1, 1].parent || error("The base rings of the two codes must be the same")
        new(A, B, repr)
    end

    function LPCode(c₁::LiftedCode, c₂::LiftedCode, repr::Function)
        c₁.A[1, 1].parent == c₂.A[1, 1].parent || error("The base rings of the two codes must be the same")
        new(c₁.A, c₂.A, repr)
    end
end

function LPCode(A::PermGroupRingMatrix, B::PermGroupRingMatrix)
    A[1, 1].parent == B[1, 1].parent || error("The base rings of the two codes must be the same")
    LPCode(A, B, permutation_repr)
end

function LPCode(c₁::LiftedCode, c₂::LiftedCode)
    c₁.A[1, 1].parent == c₂.A[1, 1].parent || error("The base rings of the two codes must be the same")
    LPCode(c₁, c₂, c₁.repr) # use the same `repr` as the first code
end

function LPCode(perm_array1::Matrix{<:Perm}, perm_array2::Matrix{<:Perm})
    LPCode(LiftedCode(perm_array1), LiftedCode(perm_array2))
end

function LPCode(shift_array1::Matrix{Int}, shift_array2::Matrix{Int}, l::Int)
    LPCode(LiftedCode(shift_array1, l), LiftedCode(shift_array2, l))
end

iscss(::Type{LPCode}) = true

function parity_checks_xz(c::LPCode)
    hx, hz = hgp(c.A, c.B')
    hx, hz = lift(c.repr, hx), lift(c.repr, hz)
    return hx, hz
end

parity_checks_x(c::LPCode) = parity_checks_xz(c)[1]

parity_checks_z(c::LPCode) = parity_checks_xz(c)[2]

parity_checks(c::LPCode) = parity_checks(CSS(parity_checks_xz(c)...))

code_n(c::LPCode) = size(c.repr(parent(c.A[1, 1])(0)), 2) * (size(c.A, 2) * size(c.B, 2) + size(c.A, 1) * size(c.B, 1))

code_s(c::LPCode) = size(c.repr(parent(c.A[1, 1])(0)), 1) * (size(c.A, 1) * size(c.B, 2) + size(c.A, 2) * size(c.B, 1))

"""
Two-block group algebra (2GBA) codes.
"""
function two_block_group_algebra_codes(a::PermGroupRingElem, b::PermGroupRingElem)
    A = reshape([a], (1, 1))
    B = reshape([b], (1, 1))
    LPCode(A, B)
end

"""
Generalized bicycle codes.
"""
function generalized_bicycle_codes(a_shifts::Array{Int}, b_shifts::Array{Int}, l::Int)
    R = PermutationGroupRing(GF(2), l)
    a = sum(R(cyclic_permutation(n, l)) for n in a_shifts)
    b = sum(R(cyclic_permutation(n, l)) for n in b_shifts)
    two_block_group_algebra_codes(a, b)
end

"""
Bicycle codes.
"""
function bicycle_codes(a_shifts::Array{Int}, l::Int)
    R = PermutationGroupRing(GF(2), l)
    a = sum(R(cyclic_permutation(n, l)) for n in a_shifts)
    two_block_group_algebra_codes(a, a')
end
