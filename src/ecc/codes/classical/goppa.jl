"""The family of Goppa codes, as discovered by Denisovich Goppa, in his 1970 paper [goppa1970new](@cite).

An irreducible binary Goppa code is characterized by two key elements:

- Goppa Polynomial `(g(x))`: A polynomial of degree `t` defined over the finite field `GF(2ᵐ)`. This polynomial must have no repeated roots. In other words, Let g(x) be a monic polynomial over the Fₓᵐ and let L={γ₀, ..., γₙ₋₁} be a set of n elements of Fₓᵐ such that g(γᵢ) ≠ 0 for 0 ≤ i < n.

- Support List `(L)`: A list containing `n` distinct elements from the finite field `GF(2ᵐ)`. These elements must not be roots of the Goppa polynomial g(x). 

The sequence L plays a crucial role in defining the structure of an irreducible binary Goppa code. The steps involved in forming L are as follows:

1. Field Selection: A finite field, `GF(2ᵐ)`, is chosen based on the extension degree `m`.
2. Goppa Polynomial: The input polynomial g(x) is used, which should be irreducible and have no repeated roots within `GF(2ᵐ)`.
3. Elements Selection for `L`: The sequence `L` consists of `n` distinct elements chosen from `GF(2ᵐ)`. Importantly, no element in L can be a root of g(x).

The parity-check matrix `(H)` of an irreducible binary Goppa code can be broken down into product of three simpler matrices:

- Vandermonde Matrix `(V)`:
  - Each row of `V` is constructed using the elements from the support list `(L)` raised to different powers.

```
1                1                  1                ...        1
(L₁)¹            (L₂)¹               (L₃)¹            ...        (Lₙ)¹ 
(L₁)²            (L₂)²               (L₃)²            ...        (Lₙ)² 
  .                .                  .                          .
  .                .                  .                          . 
  .                .                  .                          .
(L₁)ᵗ⁻¹           (L₂)ᵗ⁻¹              (L₃)ᵗ⁻¹           ...        (Lₙ)ᵗ⁻¹
```

- Diagonal Matrix `(D)`:
  - Each diagonal element in `D` is related to the roots of the Goppa polynomial `(g(x))`.

```   
1/g(L₁)		0			0			...			0
0			1/g(L₂)		0			...			0
0			0			1/g(L₃)		...			0
.			.			.			... 			.
.			.			.			... 			.
.			.			.			...			.
0			0			0			...			1/g(Lₙ)
```

- X Matrix `(X)`:
  - Each element in `X` is related to the coefficients of the Goppa polynomial `(g(x))`

```
gₜ			0			0			...			0
gₜ₋₁			gₜ			0			...			0 
gₜ₋₂			gₜ₋₁			gₜ			...			0 
.			.			.			...			.
.			.			.			...			.
.			.			.			...			.
g₁			g₂			g₃			...			gₜ
```

The matrix `VD` can be viewed as parity check matrix for `Γ(L, g)` over `GF(2ᵐ)`. The matrix is as follows:

```
1·g(L₁)⁻¹			1·g(L₂)⁻¹			1·g(L₃)⁻¹			...			1·g(Lₙ)⁻¹
(L₁)¹·g(L₁)⁻¹		(L₂)¹·g(L₂)⁻¹		(L₃)¹·g(L₃)⁻¹		...			(Lₙ)¹·g(Lₙ)⁻¹ 
(L₁)²·g(L₁)⁻¹		(L₂)²·g(L₂)⁻¹		(L₃)²·g(L₃))⁻¹		...			(Lₙ)²·g(Lₙ)⁻¹ 
	.				.				.			...				.
	.				.				.			...				. 
	.				.				.			...				.
(L₁)ᵗ⁻¹·g(L₁)⁻¹		(L₂)ᵗ⁻¹·g(L₂)⁻¹		(L₃)ᵗ⁻¹·g(L₃)⁻¹		...			(Lₙ)ᵗ⁻¹·g(Lₙ)⁻¹
```

Note: The entries of matrix over field elements are in GF(2ᵐ). Each element in GF(2ᵐ) can be represented by a `m`-tuple/binary column vector of length `m` over GF(2). If each entry of `H` is replaced by its corresponding `m`-tuple/binary column vector of length `m` over GF(2), we obtain a binary parity check matrix for the goppa code.

You might be interested in consulting [berlekamp1973goppa](@cite), [mceliece1978public](@cite), [patterson1975algebraic](@cite), [sugiyama1975method](@cite), [van1988classical](@cite), [bernstein2008attacking](@cite), [wirtz1988parameters](@cite) and [singh2019code](@cite) an as well. 

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/gappa).

Nemo Note: In Nemo, taking a random monic poly of degree `n`, this poly is irreducible with probability `1/n`. One in `n` monic polynomials is, on average irreducible. To increase probability of success of getting irreducible polynomial, use more iterations.
"""
struct Goppa <: ClassicalCode
    m::Int 
    t::Int 

    function Goppa(m, t)
        if m < 3 || t < 2 || t >= 2 ^ (m - 1)
            throw(ArgumentError("Invalid parameters: 'm' and 't' must be positive. Additionally, 'm' is >= 3 and t >= 2 to obtain a valid code and to tractable."))
        end
        new(m, t)
    end
end

function generator_polynomial(ga::Goppa)
    GF2ʳ, o = finite_field(2, ga.m, "o")
    k = GF(2, ga.m)
    po, b = polynomial_ring(k)
    gx = FqPolyRingElem
    for i in 1:500
        gx = rand(po, 1:ga.t) + b ^ ga.t
        if is_irreducible(gx) == true
            return gx  
        end
    end
end

function parity_checks(ga::Goppa)
    n = 2 ^ ga.m
    GF2ʳ, o = finite_field(2, ga.m, "o")
    gx = generator_polynomial(ga)
    L = FqFieldElem[]
    i = 0 
    while length(L) != n
        if evaluate(gx, o ^ i) != 0
            push!(L, o ^ i)
        end
        i += 1
    end
    HField = Matrix{FqFieldElem}(undef, ga.t, n)
    for i in 1:ga.t
        for j in 1:n
            HField[i, j] = L[j] ^ (i - 1) *  inv(evaluate(gx, L[j]))
        end
    end
    H = Matrix{Bool}(undef, ga.m * ga.t, n)
    for i in 1:ga.t 
        row_start = (i - 1) * ga.m + 1
        row_end = row_start + ga.m - 1
        for j in 1:n
            t_tuple = Bool[]
            for k in 0:ga.m - 1
                push!(t_tuple, !is_zero(coeff(HField[i, j], k)))
            end 
            H[row_start:row_end, j] .=  vec(t_tuple')
        end
    end 
    return H
end
