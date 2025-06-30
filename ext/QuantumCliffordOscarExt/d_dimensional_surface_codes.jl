# Oscar/AA currently does not define this and it throws error, so I have defined the fix here and submitted PR to AA.
rank(M::Oscar.Generic.DirectSumModule{T}) where T = sum(rank(summand) for summand in summands(M))

# convert from oscar mat to regular mat type
matrix_to_int(m::MatElem) = [Int(lift(ZZ, matrix(m)[i,j])) for i in 1:nrows(matrix(m)), j in 1:ncols(matrix(m))]

"""Construct the chain complex for the repetition code of length L."""
function _repcode_chain_complex(L::Int)
    F = GF(2)
    H = matrix(F, 1, L, ones(Int, L))
    V1 = free_module(F, L-1)
    V0 = free_module(F, L)
    ∂C = hom(V1, V0, H)
    return chain_complex([∂C])
end

"""Construct the chain complex for the dual of the repetition code of length L."""
function _dual_repcode_chain_complex(L::Int)
    F = GF(2)
    H = matrix(F, 1, L, ones(Int, L))
    V0 = free_module(F, L-1)
    V1 = free_module(F, L)
    ∂D = hom(V1, V0, transpose(H))
    return chain_complex([∂D])
end

# TODO:  Generalize this to D-dimensions
"""Construct a D-dimensional surface code using the hypergraph product of chain complexes.

# TODOs
# 2D surface code

# 3D surface code

# 4D surface code
"""
function d_dimensional_surface_codes(D::Int, L::Int)
    if D ∉ [2, 3, 4]
        error("Only dimensions 2, 3, and 4 are currently supported")
    end
    C = _repcode_chain_complex(L)
    D_chain = _dual_repcode_chain_complex(L)
    if D == 2
        dc = tensor_product(C, D_chain)
        E = total_complex(dc)
        Hz′ = matrix(map(E, 2))
        Hx = matrix(map(E, 1))
        Hx, Hz′ = matrix_to_int(Hx), matrix_to_int(Hz′)
        return Hx, Hz′
    elseif D == 3
        dc_2d = tensor_product(C, D_chain)
        E_2d = total_complex(dc_2d)
        dc_3d = tensor_product(E_2d, D_chain)
        F_3d = total_complex(dc_3d)
        Mz′ = matrix(map(F_3d, 3))
        Hz′ = matrix(map(F_3d, 2))
        Hx = matrix(map(F_3d, 1))
        Hx, Hz′, Mz′ = matrix_to_int(Hx), matrix_to_int(Hz′), matrix_to_int(Mz′)
        return Hx, Hz′, Mz′
    elseif D == 4
        dc_2d = tensor_product(C, D_chain)
        E_2d = total_complex(dc_2d)
        dc_3d = tensor_product(E_2d, D_chain)
        F_3d = total_complex(dc_3d)
        dc_4d = tensor_product(F_3d, C)
        G_4d = total_complex(dc_4d)
        Mz′ = matrix(map(G_4d, 4))
        Hz′ = matrix(map(G_4d, 3))
        Hx = matrix(map(G_4d, 2))
        Mx = matrix(map(G_4d, 1))
        Hx, Hz′, Mx, Mz′ = matrix_to_int(Hx), matrix_to_int(Hz′), matrix_to_int(Mx), matrix_to_int(Mz′)
        return Hx, Hz′, Mx, Mz′
    end
end
