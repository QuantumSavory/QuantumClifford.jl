# convert from oscar mat to regular mat type
matrix_to_int(m::MatElem) = [Int(lift(ZZ, matrix(m)[i,j])) for i in 1:nrows(matrix(m)), j in 1:ncols(matrix(m))]

function fq_to_int(δ::FqMatrix)
    m, n = size(δ)
    int_δ = Matrix{Int}(undef, m, n)
    @inbounds for j in 1:n, i in 1:m
        int_δ[i,j] = iszero(δ[i,j]) ? 0 : 1
    end
    return int_δ
end

const VectorFPGroupElem = Vector{FPGroupElem}

const VectorDirectProductGroupElem = Vector{GroupAlgebraElem{FqFieldElem, GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}}}

const FqFieldFPGroupAlgebra = GroupAlgebra{FqFieldElem, FPGroup, FPGroupElem}

const DirectProductGroupAlgebra = GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}
