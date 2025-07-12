# convert from oscar mat to regular mat type
matrix_to_int(m::MatElem) = [Int(lift(ZZ, matrix(m)[i,j])) for i in 1:nrows(matrix(m)), j in 1:ncols(matrix(m))]

const VectorFPGroupElem = Vector{FPGroupElem}

const VectorDirectProductGroupElem = Vector{GroupAlgebraElem{FqFieldElem, GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}}}

const FqFieldFPGroupAlgebra = GroupAlgebra{FqFieldElem, FPGroup, FPGroupElem}

const DirectProductGroupAlgebra = GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}
