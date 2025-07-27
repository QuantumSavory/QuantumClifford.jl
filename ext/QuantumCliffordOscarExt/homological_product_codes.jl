"""Convert a polynomial parity-check matrix to binary quasi-cyclic form."""
function quasi_cyclic_code(H::Generic.MatSpaceElem, l::Int)
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
Constructs a `D`-dimensional CSS quantum code (`D ≥ 2`) from `D` classical 
parity-check matrices via iterated *homological* products.

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> δ = matrix(GF(2), parity_matrix(RepCode(3)));

julia> c = HomologicalProductCode([δ,δ,δ]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(81, 3, 9)
```

Here is the `[[117, 9, 4]]` Homological product code construct from classical
quasi-cyclic code from `Table III` of [xu2024fastparallelizablelogicalcomputation](@cite)).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 3;

julia> H_poly = matrix(R, 2, 3, 
                      [x^2 x^2 x^2;
                       x   x^2  0]);

julia> H = quasi_cyclic_code(H_poly, l);

julia> c = HomologicalProductCode([H,transpose(H)]);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(117, 9, 4)
```

Here is the `[[625, 25, 9]]` Homological product code construct from classical
quasi-cyclic code from `Table III` of [xu2024fastparallelizablelogicalcomputation](@cite)).

```jldoctest
julia> using Oscar; using QuantumClifford; using QuantumClifford.ECC; using QECCore;

julia> R, x = polynomial_ring(GF(2), "x");

julia> l = 5;

julia> H_poly = matrix(R, 3, 4, 
                      [x^4   0 x^4  x^3;
                       0   x^3 x^3  x^4;
                       x^3 x^4   0  x^3]);

julia> H = quasi_cyclic_code(H_poly, l);

julia> c = HomologicalProductCode([H,transpose(H)]);

julia> code_n(c), code_k(c)
(625, 25)
```

"""
struct HomologicalProductCode{T<:MatElem} <: AbstractCSSCode
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

parity_matrix_x(hp::HomologicalProductCode) = boundary_maps(hp)[2]

parity_matrix_z(hp::HomologicalProductCode) = boundary_maps(hp)[1]'
