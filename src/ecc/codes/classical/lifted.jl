import Nemo: characteristic, lift, matrix_repr

struct LiftedCode{T} <: ClassicalCode where T<:NCRingElement
    A::Matrix{T}
    repr::Function

    function LiftedCode(A::Matrix{T}, repr::Function) where T<:NCRingElement
        new{T}(A, repr)
    end
end

function LiftedCode(A::Matrix{PermGroupRingElem{FqFieldElem}})
    !(characteristic(base_ring(A[1,1])) == 2) && error("The default permutation representation applies only to GF(2) group algebra")
    # TODO also handle the case where the base ring is not GF(2)
    LiftedCode(A, permutation_repr)
end

LiftedCode(A::Matrix{Bool}) = LiftedCode(A, x->x)

function permutation_repr(x::PermGroupRingElem{FqFieldElem})
    return sum([Int(lift(ZZ, x.coeffs[k])) .* Array(matrix_repr(k)) for k in keys(x.coeffs)], init=zeros(Bool, parent(x).l, parent(x).l))
end

function lift(repr::Function, mat::Matrix{T}) where T<:NCRingElement
    vcat([hcat([repr(mat[i, j]) for j in axes(mat, 2)]...) for i in axes(mat, 1)]...)
end

lift(A::Matrix{Bool}) = A

parity_checks(c::LiftedCode) = lift(c.repr, c.A)

code_n(c::LiftedCode) = size(c.A, 2) * size(c.repr(parent(c.A[1,1])(0)), 2)

function mod2rank(h::Matrix{<:Integer})
    Z2, _ = residue_ring(ZZ, 2)
    S = matrix_space(Z2, size(h)...)
    rank(S(h))
end

code_s(c::LiftedCode) = mod2rank(parity_checks(c)) # note that they are degenagate in general

code_k(c::LiftedCode) = code_n(c) - code_s(c)
