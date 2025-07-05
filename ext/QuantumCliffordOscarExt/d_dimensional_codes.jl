struct DDimensionalCode <: AbstractCSSCode
    """Dimension of the code (must be ≥ 2)"""
    D::Int
    """Size parameter determining the code family. For surface codes: number of physical qubits along each dimension."""
    L::Int
    """Code type, either `:Surface` or `:Toric`."""
    type::Symbol

    function DDimensionalCode(D::Int, L::Int, type::Symbol)
        D ≥ 2 || throw(ArgumentError("Dimension must be at least 2 (got D=$D)"))
        type ∈ (:Surface, :Toric) || throw(ArgumentError("Code type must be :Surface or :Toric (got :$type)"))
        new(D, L, type)
    end
end

iscss(::Type{DDimensionalCode}) = true

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

### Examples

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore: parity_matrix;

julia> D = 2; L = 2;

julia> c = d_dimensional_surface_codes(D, L);

julia> code = parity_matrix(c)
+ X_X_X
+ _X_XX
+ ZZ__Z
+ __ZZZ

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(5, 1, 2)
```

When `L = 4`, we get `[[25,1, 4]]` `2D` surface code from [Berthusen_2024](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore: parity_matrix;

julia> D = 2; L = 4;

julia> c = d_dimensional_surface_codes(D, L);

julia> code = parity_matrix(c)
+ X___X___________X________
+ _X___X__________XX_______
+ __X___X__________XX______
+ ___X___X__________X______
+ ____X___X__________X_____
+ _____X___X_________XX____
+ ______X___X_________XX___
+ _______X___X_________X___
+ ________X___X_________X__
+ _________X___X________XX_
+ __________X___X________XX
+ ___________X___X________X
+ ZZ______________Z________
+ _ZZ______________Z_______
+ __ZZ______________Z______
+ ____ZZ__________Z__Z_____
+ _____ZZ__________Z__Z____
+ ______ZZ__________Z__Z___
+ ________ZZ_________Z__Z__
+ _________ZZ_________Z__Z_
+ __________ZZ_________Z__Z
+ ____________ZZ________Z__
+ _____________ZZ________Z_
+ ______________ZZ________Z

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(25, 1, 4)
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

### [[L³ + 2L(L − 1)², 1, min(L, L²)]] 3D surface code

The `3D` surface code is obtained by taking the hypergraph product of a `2D` surface code
with a repetition code. Thus, we obtain a new `4`-term chain complex:

```math
\\begin{aligned}
F_3 \\xrightarrow{\\partial_3^F} F_2 \\xrightarrow{\\partial_2^F} F_1 \\xrightarrow{\\partial_1^F} F_0
\\end{aligned}
```

#### Metachecks

- **Z-type** metachecks: ``M_Z^T = \\partial_3^F``

### Examples

Here is an example of `[[12, 1, 2]]` `3D` Surface code with `L = 2` from [Berthusen_2024](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore: parity_matrix;

julia> D = 3; L = 2;

julia> c = d_dimensional_surface_codes(D, L);

julia> code = parity_matrix(c)
+ X_X_X_______
+ _X_XX_______
+ _____X_X_X__
+ ______X_XX__
+ X____X____X_
+ _X____X___X_
+ __X____X___X
+ ___X____X__X
+ ____X____XXX
+ ZZ__Z_____Z_
+ __ZZZ______Z
+ _____ZZ__ZZ_
+ _______ZZZ_Z

julia> code_n(c), code_k(c), distance(c)
(12, 1, 2)
```

### [[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D surface code

The `4D` surface code is constructed by taking the hypergraph product of a
`3D` surface code with a repetition code.  Thus, we obtain a new `5`-term chain complex:

```math
\\begin{aligned}
G_4 \\xrightarrow{\\partial_4^G} G_3 \\xrightarrow{\\partial_3^G} G_2 \\xrightarrow{\\partial_2^G} G_1 \\xrightarrow{\\partial_1^G} G_0
\\end{aligned}
```

### [[33, 1, 4]]

Here is an example of `[[33, 1, 4]]` `4D` Surface code with `L = 2` from [Berthusen_2024](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore: parity_matrix;

julia> D = 4; L = 2;

julia> c = d_dimensional_surface_codes(D, L);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(33, 1, 4)
```

#### Metachecks

Both X and Z-type metachecks available:
- ``M_Z^T = \\partial_4^G``
- ``M_X = \\partial_1^G``

!!! note
    To obtain surface codes of greater dimensionality, we alternate between `C` and `D` and then form a
    product with the chain complex representing the surface code in a dimension below.[Berthusen_2024](@cite).
"""
d_dimensional_surface_codes(D::Int, L::Int) = DDimensionalCode(D, L, :Surface)

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
d_dimensional_toric_codes(D::Int, L::Int) = DDimensionalCode(D, L, :Toric)

function parity_matrix_xz(c::DDimensionalCode)
    c.D >= 2 || throw(ArgumentError("Dimension must be at least 2 to construct a valid D-dimensional code."))
    hx, hz = c.type == :Surface ? _parity_matrix_xz_surface(c) : _parity_matrix_xz_toric(c)
    # Oscar returns pcm hx', so we transpose it to convert it back to hx. See page 11, B2 for reference.
    return hx', hz
end

function pcms(c::DDimensionalCode)
    c.D >= 2 || throw(ArgumentError("Dimension must be at least 2 to construct a valid D-dimensional code."))
    pcms = c.type == :Surface ? _pcms_surface(c) : _pcms_toric(c)
    return pcms
end

function _pcms_surface(c::DDimensionalCode)
    D, L = c.D, c.L
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
    boundary_maps = Vector{Matrix{Int}}()
    for d in 1:D
        ϕ = map(current, d)
        push!(boundary_maps, matrix_to_int(matrix(ϕ)))
    end
    return boundary_maps
end

function _parity_matrix_xz_surface(c::DDimensionalCode)
    boundary_maps = _pcms_surface(c)
    return c.D == 2 ? (boundary_maps[1], boundary_maps[2]) : (boundary_maps[2], boundary_maps[3])
end

parity_matrix_x(c::DDimensionalCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::DDimensionalCode) = parity_matrix_xz(c)[2]

function _pcms_toric(c::DDimensionalCode)
    c.D >= 2 || throw(ArgumentError("Dimension must be at least 2 to construct a valid D-dimensional toric code."))
    D, L = c.D, c.L
    # we use the full repetition code chain complex
    C = _repcode_chain_complex_full(L)
    current = C
    for _ in 2:D
        current = tensor_product(current, C)
        current = total_complex(current)
    end
    boundary_maps = Vector{Matrix{Int}}()
    for d in 1:D
        ϕ = map(current, d)
        push!(boundary_maps, matrix_to_int(matrix(ϕ)))
    end
    return boundary_maps
end

function _parity_matrix_xz_toric(c::DDimensionalCode)
    boundary_maps = _pcms_toric(c)
    return c.D == 2 ? (boundary_maps[1], boundary_maps[2]) : (boundary_maps[2], boundary_maps[3])
end

function parity_matrix(c::DDimensionalCode)
    c.D >= 2 || throw(ArgumentError("Dimension must be at least 2 to construct a valid D-dimensional code."))
    s = _parity_matrix_d_dimensional(c)
    return s
end

function _parity_matrix_d_dimensional(c::DDimensionalCode)
    Stabilizer(CSS(parity_matrix_xz(c)...))
end

function _chain_dimensions(C::ComplexOfMorphisms)
    rng = range(C)
    [dim(C[i]) for i in rng]
end

function code_n(c::DDimensionalCode)
    D, L = c.D, c.L
    𝒞 = _repcode_chain_complex(L)
    𝒟 = _dual_repcode_chain_complex(L)
    current = 𝒞
    for dim in 2:D
        next = dim % 2 == 0 ? 𝒟 : 𝒞
        current = tensor_product(current, next)
        current = total_complex(current)
    end
    # dimensions of all chain spaces
    dims = _chain_dimensions(current)
    # [Berthusen_2024](@cite) specifies different selection rules based on dimension, at least up to 4D.
    if D == 2
        # 2D: sum E₁ total complex dimensions
        return dims[2] # (see A11, page 9): L² + (L-1)²
    elseif D == 3
        # 3D: sum F₁ total complex dimensions
        # F₁ = E₀ ⊗ D₁ ⊕ E₁ ⊗ D₀
        return dims[3] # (see A21, page 10): L³ + 2L(L − 1)²
    elseif D == 4
        # 4D: sum G² total complex dimensions
        # G² = F₁ ⊗ C₁ ⊕ F₂ ⊗ C₀
        return dims[3] # (see A28, page 11): 6L⁴ − 12L³ + 10L² − 4L + 1
    else
        # General case: For odd D: take middle term, For even D: take middle three terms
        middle = div(length(dims), 2) + 1
        if D % 2 == 0
            return sum(dims[middle-1:middle+1])
        else
            return dims[middle]
        end
    end
end

# All D-dimensional surface codes of [Berthusen_2024](@cite) have exactly 1 logical qubit
code_k(c::DDimensionalCode) = 1

function distance(c::DDimensionalCode)
    D, L, code_type = c.D, c.L, c.type
    if code_type == :Surface
        if D == 2
            return L # 2D surface code: [[L² + (L-1)², 1, L]]
        elseif D == 3
            return min(L, L^2) # 3D surface code: [[L³ + 2L(L-1)², 1, min(L,L²)]]
        elseif D == 4
            return L^2 # 4D surface code: [[6L⁴-12L³+10L²-4L+1, 1, L²]]
        else
            throw(ArgumentError("Surface code distance not implemented for D=$D > 4. See Berthusen et al. (A4) for the general algorithm."))
        end
    else
        throw(ArgumentError("D-dimensional Toric code distance calculation not yet implemented. See Berthusen et al. (A4) for the general algorithm."))
    end
end
