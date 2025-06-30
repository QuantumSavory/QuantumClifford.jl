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

"""Construct a D-dimensional surface code using the hypergraph product of chain complexes.

# TODOs documentation
# 2D surface code

# 3D surface code

# 4D surface code
"""
function d_dimensional_surface_codes(D::Int, L::Int)
    D >= 2 || throw(ArgumentError("Dimension must be at least 2 to construct a valid D-dimensional surface code."))
    C = _repcode_chain_complex(L)
    D_chain = _dual_repcode_chain_complex(L)
    current = C
    sequence = [C]
    for dim in 2:D
        next_complex = if dim == 2
            D_chain
        elseif dim == 3
            D_chain
        else
            dim % 2 == 0 ? C : D_chain
        end
        current = tensor_product(current, next_complex)
        current = total_complex(current)
        push!(sequence, current)
    end
    boundary_maps = []
    for d in 1:D
        ϕ = map(current, d)
        push!(boundary_maps, matrix_to_int(matrix(ϕ)))
    end
    if D == 2
        return boundary_maps[1], boundary_maps[2] # Hx, Hz′
    elseif D == 3
        return boundary_maps[1], boundary_maps[2], boundary_maps[3] # Hx, Hz′, Mz′
    elseif D == 4
        return boundary_maps[2], boundary_maps[3], boundary_maps[1], boundary_maps[4] # Hx, Hz′, Mx, Mz′
    else
        # For D > 4, return all boundary maps
        return boundary_maps
    end
end

