"""We introduce a novel class of quantum CSS codes — *Multivariate Multicycle* codes — constructed
using new framework of Koszul complex over the multivariate polynomial quotient ring. For details on the construction,
please refer to our paper ([mian2026multivariatemulticyclecodescomplete](@cite))
 
Here is an example of `[[96, 44, 4]]` Multivariate Multicycle Code from Table II of [mian2026multivariatemulticyclecodescomplete](@cite).

These novel codes are made in QuantumClifford.jl backend of QuantumSavory.

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p, r = 2, 2, 2, 2;

julia> R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z]);

julia> I = ideal(R, [w^l - 1, x^m - 1, y^p - 1, z^r - 1]);

julia> S, _ = quo(R, I);

julia> A = S((1 + x)*(1 + y*z));

julia> B = S((1 + y)*(1 + z*w));

julia> C = S((1 + z)*(1 + w*x));

julia> D = S((1 + w)*(1 + x*y));

julia> c = MultivariateMulticycle([l, m, p, r], [A, B, C, D]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=900))
(96, 44, 4)
```

The Multivariate Multicycle code generalizes the several families of quantum-error correcting codes, namely:

- ## Bivariate bicycle codes ([bravyi2024high](@cite))

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

- ## Trivariate tricycle codes ([jacob2025singleshotdecodingfaulttolerantgates](@cite))

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

- ## Abelian Multicycle codes ([lin2025abelianmulticyclecodessingleshot](@cite))

Here is an example of `[[84, 6, 7]]` AMC code from Table I of [lin2025abelianmulticyclecodessingleshot](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p, r = 14, 1, 1, 1;

julia> R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z]);

julia> I = ideal(R, [w^l - 1, x^m - 1, y^p - 1, z^r - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + w);

julia> B = S(1 + w^2);

julia> C = S(1 + w^5);

julia> D = S(1 + w^6);

julia> c = MultivariateMulticycle([l, m, p, r], [A, B, C, D]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(84, 6, 7)
```

All the AMC codes from Table I are subfamilies of MM codes. Notably, this family of codes have weight-6 stabilizer checks.

- ## Cyclic Hypergraph product codes ([aydin2025cyclichypergraphproductcode](@cite))

Here is an example of `[[450, 32, 8]]` C2 code from Table I of [aydin2025cyclichypergraphproductcode](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l=15; m=15;

julia> R, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> I = ideal(R, [x^l-1, y^m-1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x + x^4);

julia> B = S(1 + y + y^4);

julia> c = MultivariateMulticycle([l,m], [A,B]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=900))
(450, 32, 8)
```

Here is an example of `[[240, 8, 8]]` CxR code from Table I of [aydin2025cyclichypergraphproductcode](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l=15; m=8;

julia> R, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> I = ideal(R, [x^l-1, y^m-1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x + x^4);

julia> B = S(1 + y);

julia> c = MultivariateMulticycle([l,m], [A,B]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=900))
(240, 8, 8)
```

- ## Abelian two-block group algebra codes ([lin2024quantum](@cite))

Here is an example of `[[16, 2, 4]]` abelian 2BGA code from Table II of [lin2024quantum](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l=2; m=4;

julia> R, (s, x) = polynomial_ring(GF(2), [:s, :x]);

julia> I = ideal(R, [s^l-1, x^m-1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x);

julia> B = S(1 + x + s + x^2 + s*x + s*x^3);

julia> c = MultivariateMulticycle([l,m], [A,B]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=900))
(16, 2, 4)
```

- ## Generalized bicycle codes ([pryadko2013quantum](@cite))

Here is an example of `[[30, 8, 4]]` generalized bicycle code from [pryadko2013quantum](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l=15; m=1;

julia> R, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> I = ideal(R, [x^l-1, y^m-1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x^2 + x^8);

julia> B = S(1 + x + x^4);

julia> c = MultivariateMulticycle([l,m], [A,B]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=900))
(30, 8, 4)
```

- ## Multivariate bicycle codes ([voss2024multivariatebicyclecodes](@cite))

Here is an example of `[[48, 4, 6]]` Weight-6 TB-QLDPC code from Appendix A Table 2 of [voss2024multivariatebicyclecodes](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l=4; m=6;

julia> R, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> I = ideal(R, [x^l-1, y^m-1]);

julia> S, _ = quo(R, I);

julia> z = x*y;

julia> A = S(x^3 + y^5);

julia> B = S(x + z^5 + y^5 + y^2);

julia> c = MultivariateMulticycle([l,m], [A,B]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=900))
(48, 4, 6)
```

- ## La-Cross codes ([pecorari2025high](@cite))

Here is an example of `[[98, 18, 4]]` "square" La-Cross from [pecorari2025high](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> n = 7;

julia> R, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> I = ideal(R, [x^n-1, y^n-1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x + x^3);

julia> B = S(1 + y + y^3);

julia> c = MultivariateMulticycle([n,n], [A,B]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=900))
(98, 18, 4)
```

- ## 4D Toric codes ([dennis2002topological](@cite))

Here is an example of `[[96, 6, 4]]` 4D Toric code from [dennis2002topological](@cite).

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p, r = 2, 2, 2, 2;

julia> R, (w, x, y, z) = polynomial_ring(GF(2), [:w, :x, :y, :z]);

julia> I = ideal(R, [w^l - 1, x^m - 1, y^p - 1, z^r - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + w);

julia> B = S(1 + x);

julia> C = S(1 + y);

julia> D = S(1 + z);

julia> c = MultivariateMulticycle([l, m, p, r], [A, B, C, D]);

julia> import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
(96, 6, 4)
```

Here is an example of a `[[1024, 30, 13 ≤ d ≤ 32]]` Haah's cubic code from Appendix B,
code D of [panteleev2021degenerate](@cite) on the `8 × 8 × 8` Lattice.

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> L = 8;

julia> R, (x,y,z) = polynomial_ring(GF(2), [:x,:y,:z]);

julia> I = ideal(R, [x^L-1, y^L-1, z^L-1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x + y + z);

julia> B = S(1 + x*y + x*z + y*z);

julia> c = MultivariateMulticycle([L,L,L], [A,B]);

julia> code_n(c), code_k(c)
(1024, 30)
```

See also: [`TrivariateTricycle`](@ref), [`BivariateBicycleViaPoly`](@ref)
"""
struct MultivariateMulticycle <: AbstractCSSCode
    orders::Vector{Int}
    polynomials::Vector{MPolyQuoRingElem{FqMPolyRingElem}}
    function MultivariateMulticycle(orders::Vector{Int}, polys::Vector{<:MPolyQuoRingElem})
        all(x->x>0, orders) || throw(ArgumentError("All orders must be positive"))
        length(polys) ≥ 2 || throw(ArgumentError("Need at least 2 variables to define a CSS code"))
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
    t = length(code.polynomials)
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
    t = length(code.polynomials)
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
    t = length(code.polynomials)
    t ≥ 4 || throw(ArgumentError("X-metachecks require t ≥ 4 variables"))
    qd = t÷2
    return transpose(maps[qd-1])
end

"""For t ≥ 3, provides metachecks for Z-stabilizers enabling single-shot decoding."""
function metacheck_matrix_z(code::MultivariateMulticycle)
    maps = boundary_maps(code)
    t = length(code.polynomials)
    t ≥ 3 || throw(ArgumentError("Z-metachecks require t ≥ 3 variables"))
    qd = t÷2
    return maps[qd+2]
end

hasmetachecks(c::MultivariateMulticycle) = length(c.polynomials) >= 4 ? (true, true) : length(c.polynomials) == 3 ? (false, true) : (false, false)
