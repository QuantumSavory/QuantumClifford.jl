"""We introduce a novel class of quantum CSS codes — *Multivariate Multicycle* codes — constructed
from multivariate polynomial quotient ring formalism over ``\\mathbb{F}_2``. 

Our discovery establishes that the boundary maps of these codes are governed by the combinatorial structure
of *Koszul* complexes. According to [eisenbud2013commutative](@cite) "Let ``a_1, \\dots, a_n`` be elements of 
``R``. Then the Koszul complex ``\\mathrm{Kosz}(\\mathbf{a})`` is *isomorphic* to the total complex of the
tensor product ``(R \\xrightarrow{a_1} R) \\otimes (R \\xrightarrow{a_2} R) \\otimes \\cdots \\otimes (R \\xrightarrow{a_n} R)``. 
For more details, see Section 17.3.

We note that the work that introduced Trivariate tricycle codes in [jacob2025singleshotdecodingfaulttolerantgates](@cite)
utilize length-1 chain complexes along with the structure of the boundary maps for the tensor-product complex of
three length-1 chain complexes that was provided in [breuckmann2024cupsgatesicohomology](@cite). See  5.3.2 Product of Λ ≥ 3 group algebra codes page 23
for more details.

Specifically, for a code defined by *t* polynomial relations, we show that the *k-th* boundary map is obtained by
taking the **Koszul matrix** in degree *k* and replacing each variable entry with the corresponding circulant matrix
derived from the code's defining relations. The **Koszul complex** provides the framework for the boundary map construction,
ensuring the commutativity properties essential for the code construction. This correspondence reveals that multivariate
multicycle codes can be constructed using the framework of **Koszul complexes**.

This family of codes generalizes the bivariate bicycle, trivariate tricycle ([`TrivariateTricycle`](@ref)), and
tetravariate tetracycle codes and it enables full single shot decoding in both X and Z directions, a capability that the
[`TrivariateTricycle`](@ref) lacks.

# Special Cases 

## t = 2: Bivariate bicycle codes ([bravyi2024high](@cite))

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l=21; m=18;

julia> R, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> I = ideal(R, [x^l-1, y^m-1]);

julia> S, _ = quo(R, I);

julia> A = S(x^3 + y^10 + y^17);

julia> B = S(y^5 + x^3  + x^19);

julia> c = MultivariateMulticycle([l,m], [A,B]);

julia> code_n(c), code_k(c)
(756, 16)
```

## t = 3: Trivariate tricycle codes ([jacob2025singleshotdecodingfaulttolerantgates](@cite))

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p = 6, 6, 4;

julia> R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]);

julia> I = ideal(R, [x^l - 1, y^m - 1, z^p - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x*y*z^3 + x^3*y^4*z^2);

julia> B = S(1 + x^3*y*z^2 + x^3*y^2*z^3);

julia> C = S(1 + x^4*y^3*z^3 + x^5*z^2);

julia> c = MultivariateMulticycle([l,m, p], [A, B, C]);

julia> code_n(c), code_k(c)
(432, 12)
```

## t = 4: Multivariate Multicycle Codes

These novel codes are made in QuantumClifford.jl backend of QuantumSavory.

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p, r = 4, 3, 3, 3;

julia> R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z]);

julia> I = ideal(R, [w^l - 1, x^m - 1, y^p - 1, z^r - 1]);

julia> S, _ = quo(R, I);

julia> A = S((1 + x^2 )*(1 + w*x*y*z^2));

julia> B = S((1 + x^2)*(1 + w*x^3*y^2*z));

julia> C = S(1 + w^2*x^2*y^2*z^2);

julia> D = S(1 + w^3*x^3*y^3*z^3);

julia> c = MultivariateMulticycle([l, m, p, r], [A, B, C, D]);

julia> code_n(c), code_k(c)
(648, 18)
```

See also: [`TrivariateTricycle`](@ref), [`BivariateBicycleViaPoly`](@ref)  
"""
struct MultivariateMulticycle <: AbstractCSSCode
    orders::Vector{Int}
    polynomials::Vector{MPolyQuoRingElem{FqMPolyRingElem}}
    function MultivariateMulticycle(orders::Vector{Int}, polys::Vector{<:MPolyQuoRingElem})
        length(orders) == length(polys) || throw(ArgumentError("Mismatch orders/polys"))
        all(x->x>0, orders) || throw(ArgumentError("All orders must be positive"))
        length(orders) ≥ 2 || throw(ArgumentError("Need at least 2 variables to define a CSS code"))
        new(orders, collect(polys))
    end
end

_gf2_to_int(M) = [iszero(M[i,j]) ? 0 : 1 for i in 1:size(M,1), j in 1:size(M,2)]

function _polynomial_to_circulant_matrix(f::MPolyQuoRingElem, orders::Vector{Int})
    t = length(orders)
    N = prod(orders)
    M = zero_matrix(GF(2), N, N)
    f_lift = lift(f)
    for col_idx in 0:(N-1)
        tmp = col_idx
        idxs = Vector{Int}(undef, t)
        for k in t:-1:1
            idxs[k] = tmp % orders[k]
            tmp ÷= orders[k]
        end
        for term in terms(f_lift)
            c = coeff(term, 1)
            iszero(c) && continue
            exps = [degree(term, k) for k in 1:t]
            row_idx = 0
            for k in 1:t
                row_idx = row_idx*orders[k]+mod(idxs[k]+exps[k], orders[k])
            end
            M[row_idx+1, col_idx+1] += c
        end
    end
    return M
end

function boundary_maps(code::MultivariateMulticycle)
    t = length(code.orders)
    N = prod(code.orders)
    @debug "Constructing boundary maps for MultivariateMulticycle"
    @debug "Number of variables t = $t"
    for k in 0:t
        dim = binomial(t, k)*N
        @debug "C_$k dimension: $dim == binom($t, $k)*$N"
    end
    circs = [_gf2_to_int(_polynomial_to_circulant_matrix(p, code.orders)) for p in code.polynomials]
    # From Wikipedia: A set of matrices A_1, ... , A_k is said to commute if they commute pairwise,
    # meaning that every pair of matrices in the set commutes. https://en.wikipedia.org/wiki/Commuting_matrices.
    # Circulant matrices commute. They form a commutative ring since the sum of two circulant matrices is circulant.
    for i in 1:t
        for j in i+1:t
            A, B = circs[i], circs[j]
            @assert mod.(A*B, 2) == mod.(B*A, 2)
        end
    end
    maps = Vector{Matrix{Int}}(undef, t)
    R, x = polynomial_ring(GF(2), ["x$i" for i in 1:t])
    K = koszul_complex(x)
    @debug "Koszul complex: $K"
    for k in 1:t
        boundary_map = map(K, k)
        KoszulMatrix = matrix(boundary_map)
        nr, nc = size(KoszulMatrix)
        for i in 1:nr, j in 1:nc
            e = KoszulMatrix[i,j]
            @assert iszero(e) || e in x "Unexpected entry in Koszul matrix at ($i,$j): $e"
        end
        M = zeros(Int, nr*N, nc*N)
        for i in 1:nr, j in 1:nc
            e = KoszulMatrix[i,j]
            !iszero(e) || continue
            found = false
            for idx in 1:t
                if e == x[idx]
                    row_range = ((i-1)*N+1):(i*N)
                    col_range = ((j-1)*N+1):(j*N)
                    M[row_range, col_range] .= circs[idx]
                    found = true
                    break
                end
            end
        end
        maps[k] = M
        @debug "Differential for degree $k: $(size(M))"
    end
    for k in 2:t
        if !isempty(maps[k]) && !isempty(maps[k-1])
            @debug "Checking exactness condition for degree $k: ∂$(k-1)∘∂$k"
            prod_check = mod.(maps[k]*maps[k-1], 2)
            @assert iszero(prod_check) "Exactness failed: ∂$(k-1)∘∂$k ≠ 0 mod 2"
        end
    end
    return maps
end

function parity_matrix_xz(code::MultivariateMulticycle)
    maps = boundary_maps(code)
    t = length(code.orders)
    if t == 2
        hx = transpose(maps[1])
        hz = maps[2]
    else
        qd = t ÷ 2
        hx = transpose(maps[qd])
        hz = maps[qd+1]
    end
    return hx, hz
end

parity_matrix_x(c::MultivariateMulticycle) = parity_matrix_xz(c)[1]

parity_matrix_z(c::MultivariateMulticycle) = parity_matrix_xz(c)[2]

"""For t ≥ 4, provides metachecks for X-stabilizers enabling single-shot decoding."""
function metacheck_matrix_x(code::MultivariateMulticycle)
    maps = boundary_maps(code)
    t = length(code.orders)
    t ≥ 4 || throw(ArgumentError("X-metachecks require t ≥ 4 variables"))
    qd = t÷2
    return transpose(maps[qd-1])
end

"""For t ≥ 3, provides metachecks for Z-stabilizers enabling single-shot decoding."""
function metacheck_matrix_z(code::MultivariateMulticycle)
    maps = boundary_maps(code)
    t = length(code.orders)
    t ≥ 3 || throw(ArgumentError("Z-metachecks require t ≥ 3 variables"))
    qd = t÷2
    return maps[qd+2]
end
