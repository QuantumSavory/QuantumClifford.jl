"""
    $TYPEDEF

Constructs a D-dimensional CSS quantum code (D ≥ 2) from D classical parity-check
matrices via iterated *homological* products.

Several interpretations of the homological product exist. For example,
[bravyi2013homologicalproductcodes](@cite) employ a simplified version known
as the *single-sector* homological product. In contrast, the `HomologicalProductCode`
adopts a more conventional definition, which [bravyi2013homologicalproductcodes](@cite)
would refer to as the *multi-sector* homological product.

The term "homological product codes" can broadly encompass various constructions
involving the product of quantum codes ([bravyi2013homologicalproductcodes](@cite),
[Campbell_2019](@cite)). However, `HomologicalProductCode` focuses specifically on a
particular subset—namely, the product of classical codes, which can also be described
as length-`1` chain complexes (sometimes called high-dimensional hypergraph product codes
[Zeng_2019](@cite)).

# Product Complex

Given ``D`` chain complexes ``\\{\\mathcal{B}^i\\}_{i\\in[D]}``, where
``\\mathcal{B}^i = \\{\\{B^i_{x_i}\\}_{x_i}, \\{\\partial^i_{x_i}\\}_{x_i}\\}``,
the ``D``-dimensional product complex is defined as:

```math
\\begin{aligned}
\\mathcal{D} = \\{\\{D_{\\vec{x}}\\}_{\\vec{x}=(x_1,\\dots,x_D)^T \\in \\mathbb{Z}^D}, \\{\\partial^i_{\\vec{x}}: D_{\\vec{x}} \\to D_{\\vec{x}-\\vec{e}_i}\\}_{i\\in[D],\\vec{x}\\in\\mathbb{Z}^D}\\} := \\text{Prod}(\\{\\mathcal{B}^i\\}_{i\\in[D]})
\\end{aligned}
```

the tensor product of these chain complexes, where:

```math
\\begin{aligned}
D_{\\vec{x}} := \\bigotimes_{i=1}^D B^i_{x_i}, \\\\
\\partial^i_{\\vec{x}} := \\bigotimes_{j=1}^D (\\partial^j_{x_j})^{\\delta_{i,j}}
\\end{aligned}
```

where ``\\delta_{i,j}`` is the Kronecker delta function and ``(\\partial^j_{x_j})^0``
is defined as the identity map [xu2024fastparallelizablelogicalcomputation](@cite).

A product complex is a high-dimensional generalization of the chain complex.

# Total Complex

The ``\\mathcal{D}``-product complex is constructed from ``D`` base chain complexes
with vector spaces ``\\{D_{\\vec{x}}\\}_{\\vec{x}}`` and boundary maps
``\\{\\partial^i_{\\vec{x}}\\}_{i\\in[D],\\vec{x}}``. It's total chain complex
``\\mathcal{T} = \\{\\{T_k\\}_k, \\{\\delta_k\\}_k\\} := \\text{Tot}(\\mathcal{D})``
as follows:

```math
\\begin{aligned}
T_k := \\bigoplus_{|\\vec{x}|=k} D_{\\vec{x}}
\\end{aligned}
```

and the boundary maps:

```math
\\begin{aligned}
\\delta_k\\left(\\bigoplus_{|\\vec{x}|=k} a_{\\vec{x}}\\right) = \\sum_{|\\vec{x}|=k} \\left(\\bigoplus_{|\\vec{y}|=k-1} \\partial_{\\vec{y},\\vec{x}} a_{\\vec{x}}\\right)
\\end{aligned}
```

for any ``a_{\\vec{x}} \\in D_{\\vec{x}}``, and

```math
\\begin{aligned}
\\partial_{\\vec{y},\\vec{x}} := 
\\begin{cases} 
\\partial^i_{\\vec{x}} & \\text{if } \\vec{x} - \\vec{y} = \\vec{e}_i \\text{ for some } i \\in [D], \\\\
0 & \\text{otherwise.}
\\end{cases}
\\end{aligned}
```

!!! note
    The total complex is obtained by projecting the ``D``-dimensional complex along
    the "diagonal" direction. Once a total chain complex is derived from a product complex
    (with length greater than 2), a quantum code can be defined from a length-2 
    subcomplex [xu2024fastparallelizablelogicalcomputation](@cite).

[xu2024fastparallelizablelogicalcomputation](@cite) focuses on product complexes with
length-1 base complexes ``\\{\\mathcal{C}^i\\}_{i\\in[D]}}`` (classical codes). In this
case, the total complex ``\\mathcal{T} = \\text{Tot}(\\text{Prod}(\\{\\mathcal{C}^i\\}_{i\\in[D]}))``
has length ``D``:

```math
\\begin{aligned}
T_D \\xrightarrow{\\delta_D} T_{D-1} \\xrightarrow{\\delta_{D-1}} \\cdots \\xrightarrow{\\delta_1} T_0
\\end{aligned}
```

[xu2024fastparallelizablelogicalcomputation](@cite) considers dimensions ``D = 2, 3, 4``.
- For ``D = 2``, the standard hypergraph product code is obtained, with planar surface codes as a special case [xu2024fastparallelizablelogicalcomputation](@cite).
- For ``D = 3`` and ``D = 4``, the construction yields 3D and 4D homological product codes, with 3D surface/toric codes serving as specific instances [xu2024fastparallelizablelogicalcomputation](@cite).

# Examples

Here is a 3D Homological product code from [Quintavalle_2021](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> δ = matrix(GF(2), parity_matrix(RepCode(3)));

julia> c = HomologicalProductCode([δ,δ,δ]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(81, 3, 3)
```

Here is the `[[117, 9, 4]]` Homological product code construct from classical
quasi-cyclic code from Table III of [xu2024fastparallelizablelogicalcomputation](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 3;

julia> H = matrix(R, 2, 3, [x^2 x^2 x^2;
                            x   x^2  0]);

julia> c = HomologicalProductCode([H,transpose(H)], l);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(117, 9, 4)
```

Here is a Homological product of `(3,4)`-classical LDPC codes.

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore; using Random; using SparseArrays

julia> μ = 2; wc = 3; wr = 4;

julia> rng = MersenneTwister(42);

julia> c = random_Gallager_ldpc(rng, μ, wc, wr);

julia> H = matrix(GF(2), c.H);

julia> c = HomologicalProductCode([H,transpose(H)]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(100, 20, 2)
```

### Fields
    $TYPEDFIELDS
"""
struct HomologicalProductCode{T<:MatElem} <: AbstractCSSCode
    """A length-`D` vector of parity-check matrices of classical error-correcting codes."""
    boundary_maps::Vector{T}
    function HomologicalProductCode(boundary_maps::Vector{T}) where {T <: MatElem}
        length(boundary_maps) >= 2 || throw(ArgumentError("At least `2` boundary maps must be provided."))
        all(base_ring(δ) == base_ring(boundary_maps[1]) for δ in boundary_maps) || throw(ArgumentError("All boundary maps must have the same base ring."))
        new{T}(boundary_maps)
    end
end

"""Convert a polynomial parity-check matrix or generator check matrix to binary quasi-cyclic form."""
function _circulant_matrix_from_quasi_cyclic_polynomial_matrix(H::MatSpaceElem, l::Int)
    F = GF(2)
    r, n = size(H)
    H_bin = zero_matrix(F, r*l, n*l)
    for i in 1:r, j in 1:n
        cfs = zeros(F, l)
        for (k, c) in enumerate(coefficients(H[i,j]))
            cfs[k] = c
        end
        circ = hcat([circshift(cfs, s) for s in 0:l-1]...)
        rows = (1:l) .+ (i-1)*l
        cols = (1:l) .+ (j-1)*l
        H_bin[rows, cols] = matrix(F, circ)
    end
    return H_bin
end

function HomologicalProductCode(boundary_maps::Vector{MatSpaceElem{FqPolyRingElem}}, l::Int)
    δₛ = [_circulant_matrix_from_quasi_cyclic_polynomial_matrix(δ, l) for δ in boundary_maps]
    HomologicalProductCode(δₛ)
end

HomologicalProductCode(boundary_maps::Vector{MatSpaceElem{FqPolyRingElem}}; l::Int) = HomologicalProductCode(boundary_maps, l)

function boundary_maps(hp::HomologicalProductCode)
    length(hp.boundary_maps) < 2 && throw(ArgumentError("`HomologicalProductCode` requires at least `2` boundary maps."))
    typeof(hp.boundary_maps) == Vector{MatSpaceElem{FqPolyRingElem}} ? [_circulant_matrix_from_quasi_cyclic_polynomial_matrix(δ, c.l) for δ in hp.boundary_maps] : hp.boundary_maps
    R = base_ring(hp.boundary_maps[1])
    C = chain_complex([hom(free_module(R, size(hp.boundary_maps[1],1)), free_module(R, size(hp.boundary_maps[1],2)), hp.boundary_maps[1])], seed=0)
    for i in 2:length(hp.boundary_maps)
        δ = hp.boundary_maps[i]
        C_next = chain_complex([hom(free_module(R, size(δ,1)), free_module(R, size(δ,2)), δ)], seed=0)
        C = tensor_product(C, C_next)
        C = total_complex(C)
    end
    δs = [matrix(map(C, d)) for d in 1:length(hp.boundary_maps)]
    return matrix_to_int.(δs)
end

function metacheck_matrix_x(c::HomologicalProductCode)
    length(c.boundary_maps) ≥ 3 || throw(ArgumentError("`X`-metachecks (`Mx`) require three classical seed codes (D=$(length(c.boundary_maps))"))
    return boundary_maps(c)[3] # Mx = δ₂
end

parity_matrix_x(hp::HomologicalProductCode) = boundary_maps(hp)[2]

parity_matrix_z(hp::HomologicalProductCode) = boundary_maps(hp)[1]'

"""
Constructs the *Double Homological Product code* from [Campbell_2019](@cite).

# 4-term Chain Complex

To construct a quantum error-correcting code with *metachecks*, we require a length-4 chain
complex. This can be built by taking the *homological product* of two length-2 chain complexes.

The length-4 chain complex, structured as:  

```math
\\begin{aligned}
\\breve{C}_{-2} \\xrightarrow{\\breve{\\delta}_{-2}} \\breve{C}_{-1} \\xrightarrow{\\breve{\\delta}_{-1}} \\breve{C}_0 \\xrightarrow{\\breve{\\delta}_0} \\breve{C}_1 \\xrightarrow{\\breve{\\delta}_1} \\breve{C}_2
\\end{aligned}
```

The homological product of two 2D chain complexes produces this length-4 complex, following the general rule:  

```math
\\begin{aligned}
\\breve{C}_m = \\bigoplus_{i - j = m} \\tilde{C}_i \\otimes \\tilde{C}_j
\\end{aligned}
```

The boundary maps are represented as block matrices and are defined as:  

```math
\\begin{align}
\\breve{\\delta}_{-2} &= \\begin{pmatrix} 
I \\otimes \\tilde{\\delta}_0^T \\\\ 
\\tilde{\\delta}_{-1} \\otimes I 
\\end{pmatrix} \\\\
\\breve{\\delta}_{-1} &= \\begin{pmatrix} 
I \\otimes \\tilde{\\delta}_{-1}^T & 0 \\\\
\\tilde{\\delta}_{-1} \\otimes I & I \\otimes \\tilde{\\delta}_0^T \\\\
0 & \\tilde{\\delta}_0 \\otimes I 
\\end{pmatrix} \\\\
\\breve{\\delta}_0 &= \\begin{pmatrix} 
\\tilde{\\delta}_{-1} \\otimes I & I \\otimes \\tilde{\\delta}_{-1}^T & 0 \\\\
0 & \\tilde{\\delta}_0 \\otimes I & I \\otimes \\tilde{\\delta}_0^T 
\\end{pmatrix} \\\\
\\breve{\\delta}_1 &= \\begin{pmatrix} 
\\tilde{\\delta}_0 \\otimes I & I \\otimes \\tilde{\\delta}_{-1}^T 
\\end{pmatrix}
\\end{align}
```

The condition ``\\breve{\\delta}_{j+1} \\breve{\\delta}_j = 0`` holds for all j,
which follows from the corresponding property of the ``\\tilde{\\delta}`` matrices.  

#### Example

Here is `[[241, 1, 9]]` double homological product code from Table I of [Campbell_2019](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> δ = [1 1 0;
            0 1 1];

julia> c = DoubleHomologicalProductCode(δ);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(241, 1, 9)
```

### Fields
    $TYPEDFIELDS
"""
struct DoubleHomologicalProductCode <: AbstractCSSCode
    """The parity-check matrix of a classical error-correcting code."""
    H::AbstractMatrix
end

function boundary_maps(c::DoubleHomologicalProductCode)
    R = GF(2)
    n = size(c.H, 2)
    m = size(c.H, 1)
    δ₀ = matrix(R, c.H)
    # see Eq. 43 of [Campbell_2019](@cite).
    # δ̃₋₁ = (I ⊗ δ₀ᵀ) ⊕ (δ₀ ⊗ I)
    δ̃₋₁ = vcat(
        kronecker_product(identity_matrix(R, n), transpose(δ₀)),
        kronecker_product(δ₀, identity_matrix(R, m))
    )
    # see Eq. 43 of [Campbell_2019](@cite).
    # δ̃₀ = [δ₀ ⊗ I | I ⊗ δ₀ᵀ]
    δ̃₀ = hcat(
        kronecker_product(δ₀, identity_matrix(R, n)),
        kronecker_product(identity_matrix(R, m), transpose(δ₀))
    )
    ñ₋₁ = n*m # Eq. 44
    ñ₀ = n^2 + m^2 # Eq. 44
    ñ₁ = n*m # Eq. 44

    # see Eq. 58 of [Campbell_2019](@cite).
    # δ̌₋₂  = (I ⊗ δ̃₀ᵀ) ⊕ (δ̃₋₁ ⊗ I)
    δ̌₋₂ = vcat(
        kronecker_product(identity_matrix(R, ñ₋₁), transpose(δ̃₀)),
        kronecker_product(δ̃₋₁, identity_matrix(R, ñ₁))
    )
    # see Eq. 59 of [Campbell_2019](@cite).
    # δ̌₋₁ = [I ⊗ δ̃₋₁ᵀ |    0
    #        δ̃₋₁ ⊗ I  | I ⊗ δ̃₀ᵀ
    #           0      | δ̃₀ ⊗ I]
    
    # I ⊗ δ̃₋₁ᵀ  0
    tb = hcat(
        kronecker_product(identity_matrix(R, ñ₁), transpose(δ̃₋₁)),
        zero_matrix(R, ñ₋₁ * ñ₋₁, ñ₁ * ñ₀)
    )

    # δ̃₋₁ ⊗ I  I ⊗ δ̃₀ᵀ
    mb = hcat(
        kronecker_product(δ̃₋₁, identity_matrix(R, ñ₀)),
        kronecker_product(identity_matrix(R, ñ₀), transpose(δ̃₀))
    )

    # 0 δ̃₀ ⊗ I
    bb = hcat(
        zero_matrix(R, ñ₁ * ñ₁, ñ₋₁ * ñ₀),
        kronecker_product(δ̃₀, identity_matrix(R, ñ₁))
    )

    δ̌₋₁ = vcat(tb, mb, bb)
    @assert iszero(δ̌₋₁*δ̌₋₂)

    # see Eq. 60 of [Campbell_2019](@cite).
    # δ̌0 = [δ̃₋₁ ⊗ I | I ⊗ δ̃₋₁ᵀ | 0
    #       0        | δ̃₀ ⊗ I   | I ⊗ δ̃₀ᵀ]
    tb = vcat(
        kronecker_product(δ̃₋₁, identity_matrix(R, ñ₋₁)),
        zero_matrix(R, ñ₀ * ñ₋₁, ñ₁ * ñ₁)
    )
    mb = vcat(
        kronecker_product(identity_matrix(R, ñ₀), transpose(δ̃₋₁)),
        kronecker_product(δ̃₀, identity_matrix(R, ñ₀))
    )
    bb = vcat(
        zero_matrix(R, ñ₀ * ñ₋₁, ñ₁ * ñ₁),
        kronecker_product(identity_matrix(R, ñ₁), transpose(δ̃₀)),
    )
    δ̌₀ = hcat(tb, mb, bb)

    # see Eq. 61 of [Campbell_2019](@cite).
    # δ̌1 = [δ̃₀ ⊗ I | I ⊗ δ̃₋₁ᵀ]
    δ̌₁ = hcat(
        kronecker_product(δ̃₀, identity_matrix(R, ñ₋₁)),
        kronecker_product(identity_matrix(R, ñ₋₁), transpose(δ̃₋₁))
    )
    @assert iszero(δ̌₁*δ̌₀)
    return fq_to_int(δ̌₋₂), fq_to_int(δ̌₋₁), fq_to_int(δ̌₀), fq_to_int(δ̌₁)
end

parity_matrix_x(c::DoubleHomologicalProductCode) = boundary_maps(c)[3] # δ̌₀

parity_matrix_z(c::DoubleHomologicalProductCode) = boundary_maps(c)[2]' # δ̌₋₁'

metacheck_matrix_x(c::DoubleHomologicalProductCode) = boundary_maps(c)[4] # δ̌₁

metacheck_matrix_z(c::DoubleHomologicalProductCode) = boundary_maps(c)[1]' # δ̌₋₂'
