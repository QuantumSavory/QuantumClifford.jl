"""Convert a polynomial parity-check matrix or generator check matrix to binary quasi-cyclic form."""
function quasi_cyclic_code(H::Oscar.Generic.MatSpaceElem, l::Int)
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

"""
    $TYPEDEF

Constructs a `D`-dimensional CSS quantum code (`D ≥ 2`) from `D` classical 
parity-check matrices via iterated *homological* products.

!!! note
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

Here is a `3D` Homological product code from [Quintavalle_2021](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> δ = matrix(GF(2), parity_matrix(RepCode(3)));

julia> c = HomologicalProductCode([δ,δ,δ]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(81, 3, 3)
```

Here is the `[[117, 9, 4]]` Homological product code construct from classical
quasi-cyclic code from `Table III` of [xu2024fastparallelizablelogicalcomputation](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 3;

julia> H_poly = matrix(R, 2, 3, [x^2 x^2 x^2;
                                 x   x^2  0]);

julia> H = quasi_cyclic_code(H_poly, l);

julia> c = HomologicalProductCode([H,transpose(H)]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(117, 9, 4)
```

Here is the `[[225, 9, 6]]` Homological product code construct from classical
quasi-cyclic code from `Table III` of [xu2024fastparallelizablelogicalcomputation](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 3;

julia> H_poly = matrix(R, 3, 4, [x^2 x^2 x^2   0;
                                 x^2   0 x^2  x^2;
                                 x^2 x^2   x  x^2]);

julia> H = quasi_cyclic_code(H_poly, l);

julia> c = HomologicalProductCode([H,transpose(H)]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(225, 9, 6)
```

Here is a Homological product of classical Reed-Muller codes.

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> r, m = 1, 3;

julia> H = matrix(GF(2), parity_matrix(ReedMuller(r,m)));

julia> c = HomologicalProductCode([H,transpose(H)]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(80, 16, 4)
```

Here is a Homological product of `(3,4)`-classical LDPC code.

julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> μ = 2; wc = 3; wr = 4;

julia> c = GallagerLDPC(μ, wc, wr);

julia> H = matrix(GF(2), parity_matrix(c));

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

function boundary_maps(hp::HomologicalProductCode)
    length(hp.boundary_maps) < 2 && throw(ArgumentError("`HomologicalProductCode` requires at least `2` boundary maps."))
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

function parity_matrix_xz(hp::HomologicalProductCode)
    hx, hz = boundary_maps(hp)[2], boundary_maps(hp)[1]'
    return hx, hz
end

function metacheck_matrix_x(c::HomologicalProductCode)
    length(c.boundary_maps) ≥ 3 || throw(ArgumentError("`X`-metachecks (`Mx`) require three classical seed codes (D=$(length(c.boundary_maps))"))
    return boundary_maps(c)[3] # Mx = δ₂
end

parity_matrix_x(hp::HomologicalProductCode) = boundary_maps(hp)[2]

parity_matrix_z(hp::HomologicalProductCode) = boundary_maps(hp)[1]'

"""
Constructs the *Double Homological Product code* from [Campbell_2019](@cite).

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

Here is `[[486, 6, 9]]` double homological product code from Table I of [Campbell_2019](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> δ = [1 1 0;
            0 1 1;
            1 0 1];

julia> c = DoubleHomologicalProductCode(δ);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(486, 6, 9)
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

function parity_matrix_xz(c::DoubleHomologicalProductCode)
    _, δ̌₋₁, δ̌₀, _ = boundary_maps(c)
    hx = δ̌₀
    hz = δ̌₋₁'
    return hx, hz
end

parity_matrix_x(c::DoubleHomologicalProductCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::DoubleHomologicalProductCode) = parity_matrix_xz(c)[2]

metacheck_matrix_x(c::DoubleHomologicalProductCode) = boundary_maps(c)[4]

metacheck_matrix_z(c::DoubleHomologicalProductCode) = boundary_maps(c)[1]'
