# Oscar/AA currently does not define this and it throws error, so I have defined the fix here and submitted PR to AA.
rank(M::Oscar.Generic.DirectSumModule{T}) where T = sum(rank(summand) for summand in summands(M))

# convert from oscar mat to regular mat type
matrix_to_int(m::MatElem) = [Int(lift(ZZ, matrix(m)[i,j])) for i in 1:nrows(matrix(m)), j in 1:ncols(matrix(m))]

const VectorFPGroupElem = Vector{FPGroupElem}

const VectorDirectProductGroupElem = Vector{GroupAlgebraElem{FqFieldElem, GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}}}

const FqFieldFPGroupAlgebra = GroupAlgebra{FqFieldElem, FPGroup, FPGroupElem}

const DirectProductGroupAlgebra = GroupAlgebra{FqFieldElem, DirectProductGroup, BasicGAPGroupElem{DirectProductGroup}}
