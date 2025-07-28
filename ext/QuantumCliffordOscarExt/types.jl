# convert from oscar mat to regular mat type
matrix_to_int(m::MatElem) = [Int(lift(ZZ, matrix(m)[i,j])) for i in 1:nrows(matrix(m)), j in 1:ncols(matrix(m))]

function fq_to_int(m::FqMatrix)
    m, n = size(m)
    H = zeros(Int, m, n)
    for i in 1:m, j in 1:n
        H[i,j] = iszero(m[i,j]) ? 0 : 1
    end
    return H
end

const VectorFPGroupElem = Vector{FPGroupElem}

const VectorDirectProductGroupElem = Vector{GroupAlgebraElem{FqFieldElem, GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}}}

const FqFieldFPGroupAlgebra = GroupAlgebra{FqFieldElem, FPGroup, FPGroupElem}

const DirectProductGroupAlgebra = GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}
