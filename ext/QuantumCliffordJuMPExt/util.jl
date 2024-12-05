function get_stab(matrix::SparseMatrixCSC{Int, Int}, logical_operator_type::Symbol)
    rows, cols = size(matrix)
    if logical_operator_type == :Z
        col_range = 1:div(cols, 2)
    else logical_operator_type == :X
        col_range = div(cols, 2) + 1:cols
    end
    submatrix = matrix[:, col_range]
    non_zero_rows = unique(submatrix.rowval)
    zero_rows = setdiff(1:rows, non_zero_rows)
    return matrix[zero_rows, :]
end

function get_lx_lz(c::Stabilizer)
    lx = stab_to_gf2(logx_ops(c))
    lz = stab_to_gf2(logz_ops(c))
    lx = SparseMatrixCSC{Int, Int}(lx)
    lz = SparseMatrixCSC{Int, Int}(lz)
    return lx, lz
end
