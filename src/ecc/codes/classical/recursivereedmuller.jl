"""
A construction of the Reed-Muller class of codes using the recursive definition.

The Plotkin `(u, u + v)` construction defines a recursive relation between generator matrices of Reed-Muller `(RM)` codes [abbe2020reed](@cite). To derive the generator matrix `G(m, r)` for `RM(r, m)`, the generator matrices of lower-order codes are utilized:
- `G(r - 1, m - 1)`: Generator matrix of `RM(r - 1, m - 1)`
- `G(r, m - 1)`: Generator matrix of `RM(r, m - 1)`

The generator matrix `G(m, r)` of `RM(m, r)` is formulated as follows in matrix notation:

```math
G(m, r) = \begin{bmatrix}
G(r, m - 1) & G(r, m - 1) \\
0 & G(r - 1, m - 1)
\end{bmatrix}
```

Here, the matrix 0 denotes an all-zero matrix with dimensions matching `G(r - 1, m - 1)`. This recursive approach facilitates the construction of higher-order Reed-Muller codes based on the generator matrices of lower-order codes.

In addition, the dimension of `RM(m - r - 1, m)` equals the dimension of the dual of `RM(r, m)`. Thus, `RM(m - r - 1, m) = RM(r, m)^⊥` shows that the [dual code](https://en.wikipedia.org/wiki/Dual_code) of `RM(r, m)` is `RM(m − r − 1, m)`, indicating the parity check matrix of `RM(r, m)` is the generator matrix for `RM(m - r - 1, m)`.

See also: `ReedMuller`

"""
struct RecursiveReedMuller <: ClassicalCode
    r::Int
    m::Int

    function RecursiveReedMuller(r, m)
        if r < 0 || r > m
            throw(ArgumentError("Invalid parameters: r must be non-negative and r ≤ m in order to valid code."))
        end
        new(r, m)
    end
end

function _recursiveReedMuller(r::Int, m::Int)
    if r == 1 && m == 1
        return Matrix{Int}([1 1; 0 1])
    elseif r == m
        return Matrix{Int}(I, 2^m, 2^m)
    elseif r == 0
        return Matrix{Int}(ones(1, 2^m))
    else
        Gᵣₘ₋₁ = _recursiveReedMuller(r, m - 1)
        Gᵣ₋₁ₘ₋₁ = _recursiveReedMuller(r - 1, m - 1)
        return vcat(hcat(Gᵣₘ₋₁, Gᵣₘ₋₁), hcat(zeros(Int, size(Gᵣ₋₁ₘ₋₁)...), Gᵣ₋₁ₘ₋₁))
    end
end

function generator(c::RecursiveReedMuller)
    return _recursiveReedMuller(c.r, c.m)
end

function parity_checks(c::RecursiveReedMuller)
    H = generator(RecursiveReedMuller(c.m - c.r - 1, c.m))
    return H
end

code_n(c::RecursiveReedMuller) = 2 ^ c.m

code_k(c::RecursiveReedMuller) = sum(binomial.(c.m, 0:c.r))

distance(c::RecursiveReedMuller) = 2 ^ (c.m - c.r)

rate(c::RecursiveReedMuller) = code_k(c::RecursiveReedMuller) / code_n(c::RecursiveReedMuller)
