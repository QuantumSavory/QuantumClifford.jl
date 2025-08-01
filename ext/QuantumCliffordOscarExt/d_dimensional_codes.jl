abstract type DDimensionalCode <: AbstractCSSCode end

"""Construct the chain complex for the repetition code of length `L`."""
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

"""Construct the chain complex for the dual of the repetition code of length `L`."""
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

"""Construct the chain complex for the repetition code of length `L`."""
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
    $TYPEDEF

Constructs the `D`-dimensional surface code using [chain complexes](https://en.wikipedia.org/wiki/Chain_complex)
and ``\\mathbb{F}_2``-[homology](https://en.wikipedia.org/wiki/Homology_(mathematics)).

## Homological Algebra Foundations of Quantum Error Correction

The theory of chain complexes over ``\\mathbb{F}_2`` provides a unified framework for
understanding error-correcting codes, where classical ``[n, k, d]`` codes correspond to
`2`-term complexes and quantum CSS codes arise naturally as `3`-term complexes satisfying
the commutativity condition ``H_Z^T H_X = 0``. This approach reveals deep connections between:

- Homological algebra and code parameters
- Boundary operators and parity check matrices
- Chain complex [exactness](https://en.wikipedia.org/wiki/Exact_sequence) and code commutativity conditions

## Chain Complex Structure

A chain complex `C` is defined by:

```math
\\begin{aligned}
C : C_n \\xrightarrow{\\partial_n} C_{n-1} \\xrightarrow{\\partial_{n-1}} \\cdots \\xrightarrow{\\partial_2} C_1 \\xrightarrow{\\partial_1} C_0
\\end{aligned}
```

with boundary operators satisfying ``\\partial_i \\circ \\partial_{i+1} = 0``. We define:

- **i-chains**: Elements of ``C_i``
- **i-cycles**: ``Z_i(C) := \\ker \\partial_i``
- **i-boundaries**: ``B_i(C) := \\mathrm{im} \\partial_{i+1}``
- **i-th homology**: ``H_i(C) := Z_i(C)/B_i(C)``

The dual complex has:

- **i-cocycles**: ``Z^i(C) := \\ker \\partial_{i+1}^T``
- **i-coboundaries**: ``B^i(C) := \\mathrm{im} \\partial_i^T``
- **i-th cohomology**: ``H^i(C) := Z^i(C)/B^i(C)``

### Classical Codes via Chain Complexes and ``\\mathbb{F_2}`` Homology

An ``[n,k,d]`` classical code corresponds to a `2`-term complex:

```math
\\begin{aligned}
0 \\rightarrow C_1 \\xrightarrow{\\partial_1 = H} C_0 \\rightarrow 0
\\end{aligned}
```

where

- ``C_1 = \\mathbb{F}_2^n`` (codeword space)
- ``C_0 = \\mathbb{F}_2^{n-k}`` (syndrome space)
- ``H`` is the parity check matrix

### Quantum CSS Codes via Chain Complexes and ``\\mathbb{F_2}`` Homology

Quantum CSS codes extend this to `3`-term complexes:

```math
\\begin{aligned}
C_2 \\xrightarrow{\\partial_2 = H_Z^T} C_1 \\xrightarrow{\\partial_1 = H_X} C_0
\\end{aligned}
```

where

- ``C_1 = \\mathbb{F}_2^n`` (physical qubits)
- ``C_2 = \\mathbb{F}_2^{m_Z}`` (`Z`-stabilizers)
- ``C_0 = \\mathbb{F}_2^{m_X}`` (`X`-stabilizers)

with the condition ``\\partial_1 \\partial_2 = H_Z^TH_X = 0`` ensuring that CSS orthogonality is satisfied.

For any chain complex, selecting two consecutive boundary operators
defines a valid CSS code. When qubits are identified with the space ``C_i``,
the code parameters are:

- number of physical qubits: ``n = \\dim C_i``
- number of logical qubits: ``k = \\dim H_i(C) = \\dim H^i(C)``
- code distance: ``d = \\min\\{\\text{wt}(v) | v \\in (H_i(C) \\cup H^i(C))\\backslash\\{0\\}\\}``

!!! note
    Quantum error-correcting codes, which are represented as `3`-term chain complexes, can be
    constructed by applying the homological or hypergraph product to two `2`-term chain complexes.

## D-dimensional Surface Code ([Berthusen_2024](@cite), [Zeng_2019](@cite))

We provide an explicit construction of the `D`-dimensional surface code
within the framework of chain complexes and homology over ``\\mathbb{F_2}``.

The quantum code is obtained by applying the homological product (or hypergraph
product) to two `2`-term chain complexes. Our construction relies on taking the
hypergraph product of these complexes.

## Double Complex

Given chain complexes `C` and `D`, we construct a [double complex](https://en.wikipedia.org/wiki/Double_complex)
derived from the tensor product of two `2`-term chain complexes:

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
\\text{Tot}(C \\boxtimes D)_i = \\bigoplus_{i=j+k} C_j \\otimes D_k = E_i
\\end{aligned}
```

with boundary maps:

```math
\\begin{aligned}
\\partial_i^E = \\bigoplus_{i=j+k} \\partial_j^v \\oplus \\partial_k^h
\\end{aligned}
```

The resulting chain complex, called the tensor product of ``C`` and ``D``, ``C ⊗ D``, enables
the construction of a CSS code when selecting any three consecutive terms in its sequence.

## Subfamilies

### [[L² + (L − 1)², 1, L]] 2D Surface Code

The `2D` surface code is constructed using the hypergraph product of two
repetition codes.Thus, we obtain a new `3`-term chain complex:

```math
\\begin{aligned}
E_2 \\xrightarrow{\\partial_2^E} E_1 \\xrightarrow{\\partial_1^E} E_0
\\end{aligned}
```

#### Examples

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC;

julia> D = 2; L = 2;

julia> c = DDimensionalSurfaceCode(D, L);

julia> code = parity_checks(c)
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
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC;

julia> D = 2; L = 4;

julia> c = DDimensionalSurfaceCode(D, L);

julia> code = parity_checks(c)
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

### [[L³ + 2L(L − 1)², 1, min(L, L²)]] 3D Surface Code

The `3D` surface code is obtained by taking the hypergraph product of a `2D` surface code
with a repetition code. Thus, we obtain a new `4`-term chain complex:

```math
\\begin{aligned}
F_3 \\xrightarrow{\\partial_3^F} F_2 \\xrightarrow{\\partial_2^F} F_1 \\xrightarrow{\\partial_1^F} F_0
\\end{aligned}
```

#### Metachecks

- **Z-type** metachecks: ``M_Z^T = \\partial_3^F``

#### Example

Here is an example of `[[12, 1, 2]]` `3D` Surface code with `L = 2` from [Berthusen_2024](@cite).

```jldoctest threeDsurface
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC;

julia> D = 3; L = 2;

julia> c = DDimensionalSurfaceCode(D, L);

julia> code = parity_checks(c)
+ XX__X_____X_
+ __XXX______X
+ _____XX__XX_
+ _______XXX_X
+ Z_Z_Z_______
+ _Z_ZZ_______
+ _____Z_Z_Z__
+ ______Z_ZZ__
+ Z____Z____Z_
+ _Z____Z___Z_
+ __Z____Z___Z
+ ___Z____Z__Z
+ ____Z____ZZZ

julia> code_n(c), code_k(c)
(12, 1)
```

!!! note
    For the `3D` surface code, there is an asymmetry between the `Z`- and `X`-bases [Berthusen_2024](@cite).
    Specifically, the `Z`-distance (``d_Z``) is `4`, whereas the `X`-distance (``d_X``) is `2`. As a result,
    the code has the parameters `[[12, 1, 2]]`.

```jldoctest threeDsurface
julia> import HiGHS; import JuMP;

julia> dz = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z))
4

julia> dx = distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X))
2
```

### [[6L⁴ − 12L³ + 10L² − 4L + 1, 1, L²]] 4D Surface Code

The `4D` surface code is constructed by taking the hypergraph product of a
`3D` surface code with a repetition code.  Thus, we obtain a new `5`-term chain complex:

```math
\\begin{aligned}
G_4 \\xrightarrow{\\partial_4^G} G_3 \\xrightarrow{\\partial_3^G} G_2 \\xrightarrow{\\partial_2^G} G_1 \\xrightarrow{\\partial_1^G} G_0
\\end{aligned}
```

#### [[33, 1, 4]]

Here is an example of `[[33, 1, 4]]` `4D` Surface code with `L = 2` from [Berthusen_2024](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC;

julia> D = 4; L = 2;

julia> c = DDimensionalSurfaceCode(D, L);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(33, 1, 4)
```

#### Metachecks

Both X and Z-type metachecks available:
- ``M_Z^T = \\partial_4^G``
- ``M_X = \\partial_1^G``

To obtain surface codes of greater dimensionality, we alternate between `C` and `D` and then form a
product with the chain complex representing the `DDimensionalSurfaceCode` [Berthusen_2024](@cite).

!!! note
    The procedure described above for the `DDimensionalSurfaceCode` can alternatively be performed
    using an `L × L` repetition code and only the chain complex `C`. In this case, the result would
    be the `DDimensionalToricCode`.

See also: [`DDimensionalToricCode`](@ref)

### Fields
    $TYPEDFIELDS
"""
struct DDimensionalSurfaceCode <: DDimensionalCode
    """Dimension of the Surface code (must be ≥ 2)."""
    D::Int
    """Size parameter determining the `D`-dimensional Surface code family, constructed via hypergraph product of 
    `(L - 1) × L` repetition code chain complexes."""
    L::Int

    function DDimensionalSurfaceCode(D::Int, L::Int)
        D ≥ 2 || throw(ArgumentError(THROW_INVALID_CODE_DIMENSION))
        new(D, L)
    end
end

function _total_complex_surface(c::DDimensionalSurfaceCode)
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
    return current
end

function _boundary_maps_surface(c::DDimensionalSurfaceCode)
    current = _total_complex_surface(c)
    boundary_maps = Vector{Matrix{Int}}()
    for d in 1:c.D
        ϕ = map(current, d)
        push!(boundary_maps, matrix_to_int(matrix(ϕ)))
    end
    return boundary_maps
end

function _parity_matrix_xz_surface(c::DDimensionalSurfaceCode)
    # 2D codes: Oscar returns pcm hx', so we transpose it to convert it back to hx. # See page 11, B2 for reference.
    # 3D and beyond codes: Oscar returns pcm hz', so we transpose it to convert it back to hz. See page 11, 12 - B3, B4, B6, and B7 for reference.
    b = _boundary_maps_surface(c)
    c.D == 2 ? (b[1]', b[2]) : (b[3], b[2]')
end

parity_matrix_xz(c::DDimensionalSurfaceCode) = _parity_matrix_xz_surface(c)

"""
    $TYPEDEF

## D-dimensional Toric Code ([Berthusen_2024](@cite), [Zeng_2019](@cite))

The `D`-dimensional toric code is obtained by taking the **iterated tensor product** of a single chain complex:

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

## Subfamilies

### 2D Toric Code

```jldoctest twoDtoric
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC;

julia> D = 2; L = 2;

julia> c = DDimensionalToricCode(D, L);

julia> code = parity_checks(c)
+ X_X_XX__
+ _X_XXX__
+ X_X___XX
+ _X_X__XX
+ ZZ__Z_Z_
+ ZZ___Z_Z
+ __ZZZ_Z_
+ __ZZ_Z_Z

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(8, 2, 2)
```

### 3D Toric Code

```jldoctest threeDtoric
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC;

julia> D = 3; L = 2;

julia> c = DDimensionalToricCode(D, L);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(24, 3, 2)
```

### 4D Toric Code

```jldoctest fourDtoric
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC;

julia> D = 4; L = 2;

julia> c = DDimensionalToricCode(D, L);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(96, 6, 4)
```

See also: [`DDimensionalSurfaceCode`](@ref)

### Fields
    $TYPEDFIELDS
"""
struct DDimensionalToricCode <: DDimensionalCode
    """Dimension of the Toric code (must be ≥ 2)."""
    D::Int
    """Size parameter determining the `D`-dimensional Toric code family, constructed via hypergraph product of 
    `L × L` repetition code chain complexes."""
    L::Int
    
    function DDimensionalToricCode(D::Int, L::Int)
        D ≥ 2 || throw(ArgumentError(THROW_INVALID_CODE_DIMENSION))
        new(D, L)
    end
end

function _total_complex_toric(c::DDimensionalToricCode)
    D, L = c.D, c.L
    # we use the full repetition code chain complex
    C = _repcode_chain_complex_full(L)
    current = C
    for _ in 2:D
        current = tensor_product(current, C)
        current = total_complex(current)
    end
    return current
end

function _boundary_maps_toric(c::DDimensionalToricCode)
    current = _total_complex_toric(c)
    boundary_maps = Vector{Matrix{Int}}()
    for d in 1:c.D
        ϕ = map(current, d)
        push!(boundary_maps, matrix_to_int(matrix(ϕ)))
    end
    return boundary_maps
end

function _parity_matrix_xz_toric(c::DDimensionalToricCode)
    # 2D codes: Oscar returns pcm hx', so we transpose it to convert it back to hx. # See page 11, B2 for reference.
    # 3D and beyond codes: Oscar returns pcm hz', so we transpose it to convert it back to hz. See page 11, 12 - B3, B4, B6, and B7 for reference.
    b = _boundary_maps_toric(c)
    c.D == 2 ? (b[1]', b[2]) : (b[3], b[2]')
end

parity_matrix_xz(c::DDimensionalToricCode) = _parity_matrix_xz_toric(c)

parity_matrix_x(c::DDimensionalCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::DDimensionalCode) = parity_matrix_xz(c)[2]

function _chain_dimensions(C::ComplexOfMorphisms)
    rng = range(C)
    [dim(C[i]) for i in rng]
end

function code_n(c::DDimensionalCode)
    current = isa(c, DDimensionalToricCode) ? _total_complex_toric(c) : _total_complex_surface(c)
    # dimensions of all chain spaces
    dims = _chain_dimensions(current)
    # [Berthusen_2024](@cite) specifies different selection rules based on dimension, at least up to 4D.
    D = c.D
    if D == 2
        # 2D: sum E₁ total complex dimensions
        return dims[2] # (see A11, page 9): For DDimensionalSurfaceCode L² + (L-1)²
    elseif D == 3
        # 3D: sum F₁ total complex dimensions
        # F₁ = E₀ ⊗ D₁ ⊕ E₁ ⊗ D₀
        return dims[2] # (see A21, page 10): For DDimensionalSurfaceCode L³ + 2L(L − 1)²
    elseif D == 4
        # 4D: sum G² total complex dimensions
        # G² = F₁ ⊗ C₁ ⊕ F₂ ⊗ C₀
        return dims[3] # (see A28, page 11): For DDimensionalSurfaceCode 6L⁴ − 12L³ + 10L² − 4L + 1
    else
        # TODO Investigate for dimension > 4? For odd D: take middle term, For even D: take middle three terms.
        middle = div(length(dims), 2) + 1
        if D % 2 == 0
            return sum(dims[middle-1:middle+1])
        else
            return dims[middle]
        end
    end
end

function code_k(c::DDimensionalCode)
    # [Berthusen_2024](@cite) uses the following left-to-right notation:
    #   D=2: E₂ → E₁ → E₀                  (qubits in E₁, compute dim H₁(E))
    #   D=3: F₃ → F₂ → F₁ → F₀       (A16) (qubits in F₁, compute dim H₁(F))
    #   D=4: G₄ → G₃ → G₂ → G₁ → G₀  (A23) (qubits in G₂, compute dim H₂(G))
    # Oscar represents the same chain complex with right-to-left notation:
    #   D=2: current_0 ← current_1 ← current_2
    #   D=3: current_0 ← current_1 ← current_2 ← current_3
    #   D=4: current_0 ← current_1 ← current_2 ← current_3 ← current_4
    current = isa(c, DDimensionalToricCode) ? _total_complex_toric(c) : _total_complex_surface(c)
    i = c.D == 2 ? 1 : 2 # For codes with D ≥ 3: Hx is boundary_maps(c)[3], boundary_maps(c)[2] is Hz'.
    # TODO Does this hold for dimension > 4?
    ∂ᵢ = map(current, i)
    ∂ᵢ₊₁ = map(current, i+1)
    ker_∂ᵢ = kernel(∂ᵢ)[1]
    im_∂ᵢ₊₁ = image(∂ᵢ₊₁)[1]
    dim_ker = dim(ker_∂ᵢ)
    dim_im = dim(im_∂ᵢ₊₁)
    # number of logical qubits == dim ker ∂ᵢ - dim im ∂ᵢ₊₁
    k = dim_ker - dim_im
    return k
end

function metacheck_matrix_x(c::DDimensionalCode)
    c.D ≥ 4 || throw(ArgumentError(THROW_INVALID_X_METACHECKS))
    return Matrix(boundary_maps(c)[4]) # Mx
end

function metacheck_matrix_z(c::DDimensionalCode)
    c.D ≥ 3 || throw(ArgumentError(THROW_INVALID_Z_METACHECKS))
    return Matrix(boundary_maps(c)[1]') # Mz
end

"""
$TYPEDEF

Returns all boundary maps of the chain complex, including both *parity check*
and *metacheck* matrices.

Here are the boundary maps of `[[12, 1, 2]]` `3D` Surface code with
`L = 2` from [Berthusen_2024](@cite).

```jldoctest boundarymaps
julia> using Oscar; using QuantumClifford; using QECCore;

julia> using QuantumClifford.ECC: DDimensionalSurfaceCode, boundary_maps, metacheck_matrix_z;

julia> D = 3; L = 2;

julia> c = DDimensionalSurfaceCode(D, L);

julia> Mz, Hz, Hx = boundary_maps(c);
```

The parity check matrices of `[[12, 1, 2]]` `3D` Surface code are

```jldoctest boundarymaps
julia> Hx
4×12 Matrix{Int64}:
 1  1  0  0  1  0  0  0  0  0  1  0
 0  0  1  1  1  0  0  0  0  0  0  1
 0  0  0  0  0  1  1  0  0  1  1  0
 0  0  0  0  0  0  0  1  1  1  0  1
```

!!! note
    For `3D` and higher-dimensional codes, `Oscar` returns `Z`-type parity check
    matrix as transpose (``H_Z^T``). We transpose it to convert it back to ``H_Z``.
    See `B3`, page `11` of [Berthusen_2024](@cite).

```jldoctest boundarymaps
julia> Hz'
9×12 adjoint(::Matrix{Int64}) with eltype Int64:
 1  0  1  0  1  0  0  0  0  0  0  0
 0  1  0  1  1  0  0  0  0  0  0  0
 0  0  0  0  0  1  0  1  0  1  0  0
 0  0  0  0  0  0  1  0  1  1  0  0
 1  0  0  0  0  1  0  0  0  0  1  0
 0  1  0  0  0  0  1  0  0  0  1  0
 0  0  1  0  0  0  0  1  0  0  0  1
 0  0  0  1  0  0  0  0  1  0  0  1
 0  0  0  0  1  0  0  0  0  1  1  1
```

Along with the `Z`- and `X`-type parity check matrices, we have a metacheck matrix
specifically for the `Z`-type checks. The classical code derived from this metacheck
matrix has a distance of `d = 2` meaning it can identify (but not correct) a single error
in the `Z`-type syndrome measurements. See page `12` of [Berthusen_2024](@cite) for details.

```jldoctest boundarymaps
julia> Mz'
2×9 adjoint(::Matrix{Int64}) with eltype Int64:
 1  0  1  0  1  0  1  0  1
 0  1  0  1  0  1  0  1  1
```

We can use `metacheck_matrix_z` directly instead of using `boundary_maps`.

```jldoctest boundarymaps
julia> metacheck_matrix_z(c)
2×9 Matrix{Int64}:
 1  0  1  0  1  0  1  0  1
 0  1  0  1  0  1  0  1  1
```

# Metachecks in CSS Codes

The parity-check matrices ``M_Z`` and ``M_X`` are called *metachecks* in CSS
codes. These matrices emerge from the constraints imposed by boundary maps, which
satisfy the condition ``\\partial_{i+1} \\partial_i = 0``. This guarantees that:

```math
\\begin{aligned}
M_Z H_Z = 0 \\quad \\text{and} \\quad M_X H_X = 0,
\\end{aligned}
```

meaning that:

- Valid ``Z``-type syndromes must be in ``\\ker M_Z``
- Valid ``X``-type syndromes must be in ``\\ker M_X``

When measurement errors occur, they distort the syndrome vector ``\\mathbf{s}``,
generating a detectable *metasyndrome*. By examining ``\\mathbf{m}``, we can identify
and correct errors in ``\\mathbf{s}`` before proceeding with standard decoding. This
technique is called *syndrome repair decoding* [Higgott_2023](@cite).
"""
boundary_maps(c::DDimensionalSurfaceCode) = _boundary_maps_surface(c)
boundary_maps(c::DDimensionalToricCode) = _boundary_maps_toric(c)
