function get_stab_hx(matrix::SparseMatrixCSC{T, Int}) where T
    rows, cols = size(matrix)
    rhs_start = div(cols, 2) + 1
    rhs_cols = matrix[:, rhs_start:cols]
    non_zero_rows_rhs = unique(rhs_cols.rowval)
    zero_rows_rhs = setdiff(1:rows, non_zero_rows_rhs)
    return matrix[zero_rows_rhs, :]
end

function get_lx_lz(c::Stabilizer)
    lx = stab_to_gf2(logicalxview(canonicalize!(MixedDestabilizer(c))))
    lz = stab_to_gf2(logicalzview(canonicalize!(MixedDestabilizer(c))))
    return lx, lz
end
