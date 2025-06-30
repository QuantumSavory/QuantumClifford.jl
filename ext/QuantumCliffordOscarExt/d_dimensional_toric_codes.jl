"""Construct the chain complex for the repetition code of length L."""
function _repcode_chain_complex_full(L::Int)
    F = GF(2)
    H = parity_matrix(RepCode(L))
    H = matrix(F, H)
    V1 = free_module(F, L)
    V0 = free_module(F, L)
    ∂C = hom(V1, V0, H)
    return chain_complex([∂C])
end

"""
The D-dimensional toric code is obtained by taking the **iterated tensor product** of a single chain complex:

```math
\\begin{aligned}
C = \\left( \\mathbb{F}_2^L \\xrightarrow{H} \\mathbb{F}_2^L \\right)
\\end{aligned}
```

where `H` is the ``L \\times L`` parity check matrix of the repetition code. The total complex
is built by taking the tensor product ``C^{\\otimes D}`` and forming the associated total complex via direct sums.

!!! note
    D-dimensional toric code construction differs from the surface code as we
    use the full ``L \\times L`` repetition code (rather than ``(L-1) \\times L``). The
    tensor products with identical complexes (rather than alternating complexes) are used.
"""
function d_dimensional_toric_codes(D::Int, L::Int)
    D >= 2 || throw(ArgumentError("Dimension must be at least 2 to construct a valid D-dimensional toric code."))
    # we use the full repetition code chain complex
    C = _repcode_chain_complex_full(L)
    current = C
    for _ in 2:D
        current = tensor_product(current, C)
        current = total_complex(current)
    end
    boundary_maps = []
    for d in 1:D
        ϕ = map(current, d)
        push!(boundary_maps, matrix_to_int(matrix(ϕ)))
    end
    if D == 2
        return boundary_maps[1], boundary_maps[2] # Hx, Hz
    elseif D == 3
        return boundary_maps[1], boundary_maps[2], boundary_maps[3] # Hx, Hz′, Mz′
    elseif D == 4
        return boundary_maps[1], boundary_maps[2], boundary_maps[4], boundary_maps[3] # Hx, Hz′, Mx, Mz′
    else
        # For D > 4, return all boundary maps
        return boundary_maps
    end
end
