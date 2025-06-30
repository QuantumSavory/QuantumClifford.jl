"""Construct the chain complex for the repetition code of length L."""
function _repcode_chain_complex(L::Int)
    F = GF(2)
    H = parity_matrix(RepCode(L))
    H = H[1:L-1, :] # (L-1) × L (A5)
    H = matrix(F, H)
    V1 = free_module(F, L-1)
    V0 = free_module(F, L)
    ∂C = hom(V1, V0, H)
    return chain_complex([∂C])
end

"""Construct the chain complex for the dual of the repetition code of length L."""
function _dual_repcode_chain_complex(L::Int)
    F = GF(2)
    H = parity_matrix(RepCode(L))
    H = H[1:L-1, :] # (L-1) × L (A5)
    H = matrix(F, H)
    V0 = free_module(F, L-1)
    V1 = free_module(F, L)
    ∂D = hom(V1, V0, transpose(H))
    return chain_complex([∂D])
end

"""

# D-dimensional Surface Code

The construction uses chain complexes and ``\\mathbb{F}_2``-homology. For a quantum
code, we take the hypergraph product of two `2`-term chain complexes.

## Double Complex

Given chain complexes `C` and `D`, we form a double complex:

```math
\\begin{aligned}
C \\boxtimes D \\quad \\text{with} \\quad \\partial_i^v = \\partial_i^C \\otimes I_{D_i} \\quad \\text{and} \\quad \\partial_i^h = I_{C_i} \\otimes \\partial_i^D
\\end{aligned}
```

## Total Complex

The total complex is derived from a double complex by taking the direct
sum of vector spaces and boundary maps that share the same dimension:

```math
\\begin{aligned}
\text{Tot}(C \\boxtimes D)_i = \\bigoplus_{i=j+k} C_j \\otimes D_k = E_i
\\end{aligned}
```

with boundary maps:

```math
\\begin{aligned}
\\partial_i^E = \\bigoplus_{i=j+k} \\partial_j^v \\oplus \\partial_k^h
\\end{aligned}
```

## Subfamilies

### [[L² + (L − 1)², 1, L]] 2D surface code

The `2D` surface code is constructed using the hypergraph product of two
repetition codes.Thus, we obtain a new `3`-term chain complex:

```math
\\begin{aligned}
E_2 \\xrightarrow{\\partial_2^E} E_1 \\xrightarrow{\\partial_1^E} E_0
\\end{aligned}
```

#### Chain Complex

The construction is as follows:

```math
\\begin{aligned}
C = \\left( C_1 \\xrightarrow{\\partial} C_0 \\right) \\quad \\text{and} \\quad D = \\left( D_1 \\xrightarrow{\\partial^T} D_0 \\right)
\\end{aligned}
```

where ``\\partial`` is the ``(L-1) \\times L`` parity check matrix:

```math
\\begin{aligned}
H = \\begin{pmatrix}
1 & 1 & & \\\\
 & 1 & \\ddots & \\\\
 & & \\ddots & 1 \\\\
 & & & 1
\\end{pmatrix}
\\end{aligned}
```

### [L³ + 2L(L − 1)², 1, min(L, L²)]] 3D surface code

The `3D` surface code is obtained by taking the hypergraph product of a `2D` surface code
with a repetition code. Thus, we obtain a new `4`-term chain complex:

```math
\\begin{aligned}
F_3 \\xrightarrow{\\partial_3^F} F_2 \\xrightarrow{\\partial_2^F} F_1 \\xrightarrow{\\partial_1^F} F_0
\\end{aligned}
```

#### Metachecks:

- **Z-type** metachecks: ``M_Z^T = \\partial_3^F``

### [[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D surface code

The `4D` surface code is constructed by taking the hypergraph product of a
`3D` surface code with a repetition code.  Thus, we obtain a new `5`-term chain complex:

```math
\\begin{aligned}
G_4 \\xrightarrow{\\partial_4^G} G_3 \\xrightarrow{\\partial_3^G} G_2 \\xrightarrow{\\partial_2^G} G_1 \\xrightarrow{\\partial_1^G} G_0
\\end{aligned}
```

#### Metachecks:

Both X and Z-type metachecks available:
- ``M_Z^T = \\partial_4^G``
- ``M_X = \\partial_1^G``

!!! note
    To obtain surface codes of greater dimensionality, we alternate between `C` and `D` and then form a
    product with the chain complex representing the surface code in a dimension below.[Berthusen_2024](@cite).
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
    # TODO: Maybe it's too soon to return the boundary maps as Int matrices, we can do Oscar operations on them.
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
