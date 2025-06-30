# Oscar/AA currently does not define this and it throws error, so I have defined the fix here.
rank(M::Generic.DirectSumModule{T}) where T = sum(rank(summand) for summand in summands(M))

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
        H_Z_T = matrix(map(E, 2))
        H_X = matrix(map(E, 1))
        return H_X, H_Z_T
    elseif D == 3
        dc_2d = tensor_product(C, D_chain)
        E_2d = total_complex(dc_2d)
        dc_3d = tensor_product(E_2d, D_chain)
        F_3d = total_complex(dc_3d)
        M_Z_T = matrix(map(F_3d, 3))
        H_Z_T = matrix(map(F_3d, 2))
        H_X = matrix(map(F_3d, 1))
        return H_X, H_Z_T, M_Z_T
    elseif D == 4
        dc_2d = tensor_product(C, D_chain)
        E_2d = total_complex(dc_2d)
        dc_3d = tensor_product(E_2d, D_chain)
        F_3d = total_complex(dc_3d)
        dc_4d = tensor_product(F_3d, C)
        G_4d = total_complex(dc_4d)
        M_Z_T = matrix(map(G_4d, 4))
        H_Z_T = matrix(map(G_4d, 3))
        H_X = matrix(map(G_4d, 2))
        M_X = matrix(map(G_4d, 1))
        n = size(H_X, 2)
        return H_X, H_Z_T, M_X, M_Z_T
    end
    # Generalize this to D-dimensions
end
