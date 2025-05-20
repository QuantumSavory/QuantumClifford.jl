function get_lx_lz(c::Stabilizer)
    lx = stab_to_gf2(logx_ops(c))
    lz = stab_to_gf2(logz_ops(c))
    lx = SparseMatrixCSC{Int, Int}(lx)
    lz = SparseMatrixCSC{Int, Int}(lz)
    return lx, lz
end
