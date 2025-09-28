
"""
We introduce a novel class of quantum CSS codes — *Multivariate Multicycle* codes — constructed
from multivariate polynomial quotient ring formalism over finite fields over GF(2). Our discovery establishes
that the boundary maps of these codes are governed by the combinatorial structure of *Koszul* complexes.
Specifically, for a code defined by *t* polynomial relations, we demonstrate that the *k-th* boundary map
is obtained by taking the **Koszul matrix** in degree *k* and replacing each variable entry with the
corresponding circulant matrix derived from the code's defining relations. The **Koszul complex** provides
the framework for the boundary map construction, ensuring the commutativity properties essential for the code
construction. This correspondence reveals that multivariate multicycle codes possess the structure of
**Koszul complexes** with circulant coefficients. This family of codes generalizes the bivariate bicycle,
trivariate tricycle ([TrivariateTricycleCode](@ref)), and tetravariate tetracycle codes and it enables
full single shot decoding in both X and Z directions.

# Special Cases 

## t = 2: Bivariate bicycle codes

```jldoctest
julia> l=21; m=18;

julia> R, (x, y) = polynomial_ring(GF(2), [:x, :y]);

julia> I = ideal(R, [x^l-1, y^m-1]);

julia> S, _ = quo(R, I);

julia> A = S(x^3 + y^10 + y^17);

julia> B = S(y^5 + x^3  + x^19);

julia> c = MultivariateMulticycleCode([l,m], [A,B])

julia> code_n(c), code_k(c)
(756, 16)
```

## t = 3: Trivariate tricycle codes

```jldoctest
julia> using Oscar; using QuantumClifford.ECC;

julia> l, m, p = 6, 6, 4;

julia> R, (x, y, z) = polynomial_ring(GF(2), [:x, :y, :z]);

julia> I = ideal(R, [x^l - 1, y^m - 1, z^p - 1]);

julia> S, _ = quo(R, I);

julia> A = S(1 + x*y*z^3 + x^3*y^4*z^2);

julia> B = S(1 + x^3*y*z^2 + x^3*y^2*z^3);

julia> C = S(1 + x^4*y^3*z^3 + x^5*z^2);

julia> c = MultivariateMulticycleCode([l,m, p], [A, B, C]);

julia> code_n(c), code_k(c)
(432, 12)
```

## t = 4: Tetravariate Tetracycle Codes

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

julia> c = MultivariateMulticycleCode([l, m, p, r], [A, B, C, D]);

julia> code_n(c), code_k(c)
(648, 18)
```
"""
struct MultivariateMulticycleCode <: AbstractCSSCode
    orders::Vector{Int}
    polynomials::Vector{MPolyQuoRingElem{FqMPolyRingElem}}
    function MultivariateMulticycleCode(orders::Vector{Int}, polys::Vector{<:MPolyQuoRingElem})
        length(orders) == length(polys) || throw(ArgumentError("Mismatch orders/polys"))
        all(x->x>0, orders) || throw(ArgumentError("All orders must be positive"))
        length(orders) ≥ 2 || throw(ArgumentError("Need at least 2 variables to define an CSS code"))
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
                row_idx = row_idx * orders[k] + mod(idxs[k] + exps[k], orders[k])
            end 
            M[row_idx+1, col_idx+1] += c
        end
    end
    return M
end

function boundary_maps(code::MultivariateMulticycleCode)
    t = length(code.orders)
    N = prod(code.orders)
    circs = [_gf2_to_int(_polynomial_to_circulant_matrix(p, code.orders)) for p in code.polynomials]
    maps = Vector{Matrix{Int}}(undef, t)
    for k in 1:t
        R, x = polynomial_ring(GF(2), ["x$i" for i in 1:t])
        KoszulMatrix = koszul_matrix(x, k)
        nr, nc = size(KoszulMatrix)
        M = zeros(Int, nr*N, nc*N)
        for i in 1:nr, j in 1:nc
            e = KoszulMatrix[i, j]
            !iszero(e) || continue
            for idx in 1:t
                e == x[idx] || continue
                r = ((i-1)*N+1):(i*N)
                c = ((j-1)*N+1):(j*N)
                M[r, c] .= circs[idx]
                break
            end
        end
        maps[k] = M
    end
    return maps
end

function parity_matrix_xz(c::MultivariateMulticycleCode; qubit_degree::Union{Nothing,Int}=nothing)
    maps = boundary_maps(c)
    t = length(c.orders)
    k = qubit_degree === nothing ? cld(t, 2) : qubit_degree
    k >= 1 || throw(ArgumentError("qubit degree must be >= 1"))
    (k+1) <= t || throw(ArgumentError("qubit degree too large for number of variables"))
    hx, hz = transpose(maps[k]), maps[k+1]
    return hx, hz
end

parity_matrix_x(c::MultivariateMulticycleCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::MultivariateMulticycleCode) = parity_matrix_xz(c)[2]
