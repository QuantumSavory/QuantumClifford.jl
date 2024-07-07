using Nemo

struct LiftedCode{T} <: ClassicalCode
    A::Union{MatrixElem{T}, Matrix{T}} where T<:NCRingElement
    repr::Function

    function LiftedCode(A::Union{MatrixElem{T}, Matrix{T}} where T<:NCRingElement , repr::Function)
        new{T}(A, repr)
    end

    function LiftedCode(A::Union{MatrixElem{PermGroupRingElem{FqFieldElem}}, Matrix{PermGroupRingElem{FqFieldElem}}})
        # permutation representation applies only to GF(2) group algebra
        @assert characteristic(A.base_ring.base_ring) == 2
        new{PermGroupRingElem{FqFieldElem}}(A, permutation_repr)
    end
end

function permutation_repr(x::PermGroupRingElem{FqFieldElem})
    return sum([Int(lift(ZZ, x.coeffs[k])) .* Array(matrix_repr(k)) for k in keys(x.coeffs)], init=zeros(Bool, parent(x).l, parent(x).l))
end

function parity_checks(c::LiftedCode)
    vcat([hcat([c.repr(c.A[i, j]) for j in 1:size(c.A, 2)]...) for i in 1:size(c.A, 1)]...)
end

code_n(c::LiftedCode) = size(parity_checks(c), 2)

code_n(c::LiftedCode{PermGroupRingElem{FqFieldElem}}) = characteristic(c.A.base_ring.base_ring) == 2 ? size(c.A, 2) * c.A.base_ring.l : size(c.A, 2)

code_s(c::LiftedCode) = rank(parity_checks(c))

code_k(c::LiftedCode) = code_n(c) - code_s(c)
