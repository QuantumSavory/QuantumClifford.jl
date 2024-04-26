"""The family of Goppa codes, as discovered by Denisovich Goppa, in his 1970 paper [goppa1970new](@cite).

An irreducible binary Goppa code is characterized by two key elements:

- Goppa Polynomial `(g(x))`: A polynomial of degree `t` defined over the finite field `GF(2^m)`. This polynomial must have no repeated roots. In other words, Let g(x) be a monic polynomial over the F_q ^ m and let L={γ_0, ..., γ_(n − 1)} be a set of n elements of F_q ^ m such that g(γ_i) ≠ 0 for 0 ≤ i < n.

- Support List `(L)`: A list containing `n` distinct elements from the finite field `GF(2^m)`. These elements must not be roots of the Goppa polynomial g(x). 

The sequence L plays a crucial role in defining the structure of an irreducible binary Goppa code. The steps involved in forming L are as follows:

1. Field Selection: A finite field, `GF(2^m)`, is chosen based on the extension degree `m`.
2. Goppa Polynomial: The input polynomial g(x) is used, which should be irreducible and have no repeated roots within `GF(2^m)`.
3. Elements Selection for `L`: The sequence `L` consists of `n` distinct elements chosen from `GF(2^m)`. Importantly, no element in L can be a root of g(x).

The parity-check matrix `(H)` of an irreducible binary Goppa code can be expressed as the product `(H = XVD)` of three simpler matrices:

- Vandermonde Matrix `(V)`:
  - Each row of `V` is constructed using the elements from the support list `(L)` raised to different powers.

1                    1                       1                    ...        1
(L_1) ^ 1            (L_2) ^ 1               (L_3) ^ 1            ...        (L_n) ^ 1 
(L_1) ^ 2            (L_2) ^ 2               (L_3) ^ 2            ...        (L_n) ^ 2 
      .                    .                       .                             .
      .                    .                       .                             . 
      .                    .                       .                             .
(L_1) ^ (t - 1)      (L_2) ^ (t - 1)         (L_3) ^ (t - 1)      ...        (L_n) ^ (t - 1)         
   
- Diagonal Matrix `(D)`:
  - Each diagonal element in `D` is related to the roots of the Goppa polynomial `(g(x))`.
      
   1/g(L_1)                0                       0              ...            0
      0                 1/g(L_2)                   0              ...            0
      0                    0                    1/g(L_3)          ...            0
      .                    .                       .                             . 
      .                    .                       .                             . 
      .                    .                       .                             .
      0                    0                       0              ...         1/g(L_n)

- X Matrix `(X)`:
  - Each element in `X` is related to the coefficients of the Goppa polynomial `(g(x))`

     g_t                   0                       0              ...            0
  g_(t - 1)               g_t                      0              ...            0 
  g_(t - 2)            g_(t - 1)                  g_t             ...            0 
      .                    .                       .              ...            .
      .                    .                       .              ...            .
      .                    .                       .              ...            .
     g_1                  g_2                     g_3             ...           g_t

You might be interested in consulting [berlekamp1973goppa](@cite), [mceliece1978public](@cite), [patterson1975algebraic](@cite), [sugiyama1975method](@cite), [van1988classical](@cite), [bernstein2008attacking](@cite), [wirtz1988parameters](@cite) and [singh2019code](@cite) an as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/gappa)

Nemo Note: In Nemo, taking a random monic poly of degree `n`, this poly is irreducible with probability `1/n`. One in `n` monic polynomials is, on average irreducible. To increase probability of success of getting irreducible polynomial, use more iterations.
"""
struct Goppa <: ClassicalCode
    n::Int 
    t::Int 

    function Goppa(n, t)
        if n < 7 || n > 500 || t < 2 || t > n
            throw(ArgumentError("Invalid parameters: 'n' and 't' must be positive. Additionally, 'n' is >= to 7 and t >= 2 and t < n to obtain a valid code and to tractable."))
        end
        new(n, t)
    end
end

function generator_polynomial(ga::Goppa)
    r = ceil(Int, log2(ga.n))
    GF2ͬ, o = finite_field(2, r, "o")
    k = GF(2, r)
    po, b = polynomial_ring(k)
    gx = FqPolyRingElem
    for i in 1:100
        gx = rand(po, 1:ga.t) + b ^ ga.t
        if is_irreducible(gx) == true
            return gx  
        end
    end
end

function parity_checks(ga::Goppa)
    r = ceil(Int, log2(ga.n))
    GF2ͬ, o = finite_field(2, r, "o")
    gx = generator_polynomial(ga)
    L = FqFieldElem[]
    i = 0 
    while length(L) != ga.n
        if evaluate(gx, o ^ i) != 0
            L = [L; evaluate(gx, o^i)]
        end
        i += 1
    end
    V = Matrix{FqFieldElem}(undef, ga.t, ga.n)
    for j in 1:ga.n
        V[1, j] = o ^ 0
    end
    for i in 2:ga.t
        for j in 1:ga.n
            V[i, j] = L[j] ^ (i - 1)
        end
    end
    M = identity_matrix(GF2ͬ, ga.n)
    for i in 1:ga.n
        for j in 1:ga.n
            M[i, j] = getindex(M, i, j)*inv(evaluate(gx, L[j]))
        end
    end
    HField = Matrix{FqFieldElem}(undef, ga.t, ga.n)
    for i in 1:ga.t
        for j in 1:ga.n
            HField[i, j] = V[i, j] * M[i, j]
        end
    end
    X = Matrix{FqFieldElem}(undef, ga.t, ga.t)
    fill!(X, o^0 - 1)
    j = 1
    for col in eachcol(X)
        col[j:ga.t] .= (coeff(gx, i) for i in ga.t - j + 1:-1:1)
        j += 1   
    end  
    HF = Matrix{FqFieldElem}(undef, ga.t, ga.n)
    HF = X * HField
    H = Matrix{Bool}(undef, r*ga.t, ga.n)
    for i in 1:ga.t 
        row_start = (i - 1) * r + 1
        row_end = row_start + r - 1
        for j in 1:ga.n
            t_tuple = Bool[]
            for k in 0:r - 1
                t_tuple = [t_tuple; !is_zero(coeff(HF[i, j], k))]
            end 
            H[row_start:row_end, j] .=  vec(t_tuple')
        end
    end 
    return H
end
