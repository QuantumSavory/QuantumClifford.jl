"""The family of Reed-Muller codes, as discovered by Muller in his 1954 paper [muller1954application](@cite) and Reed who proposed the first efficient decoding algorithm [reed1954class](@cite).

Let `m` be a positive integer and `r` a nonnegative integer with `r ≤ m`. These linear codes, denoted as `RM(r, m)`, have order `r` (where `0 ≤ r ≤ m`) and codeword length `n` of `2ᵐ`.

Two special cases of generator(RM(r, m)) exist:
    1. `generator(RM(0, m))`: This is the `0ᵗʰ`-order `RM` code, similar to the binary repetition code with length `2ᵐ`. It's characterized by a single basis vector containing all ones.
    2. `generator(RM(m, m))`: This is the `mᵗʰ`-order `RM` code. It encompasses the entire field `F(2ᵐ)`, representing all possible binary strings of length `2ᵐ`.

You might be interested in consulting [raaphorst2003reed](@cite), [abbe2020reed](@cite), and [djordjevic2021quantum](@cite) as well.

The dimension of `RM(m - r - 1, m)` equals the dimension of the dual of `RM(r, m)`. Thus, `RM(m - r - 1, m) = RM(r, m)^⊥` which indicates that the parity check matrix of `RM(r, m)` is the generator matrix for `RM(m - r - 1, m)`.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/reed_muller).
"""

abstract type ClassicalCode end

struct ReedMuller <: ClassicalCode
    r::Int
    m::Int

    function ReedMuller(r, m)
        0 ≤ r ≤ m && m - r - 1 ≥ 0 || m < 11 || throw(ArgumentError("Invalid parameters: r must be non-negative and r ≤ m. Additionally, m must be positive and < 11 in order to obtain a valid code and to remain tractable"))
        new(r, m)
    end
end

function _variablesₓᵢ_rm(m, i)
    return repeat([fill(1, 2 ^ (m - i - 1)); fill(0, 2 ^ (m - i - 1))], outer = 2 ^ i)
end

function _vmult_rm(vecs...)
    return [reduce(*, a, init=1) for a in zip(vecs...)]
end

function generator(c::ReedMuller)
    r = c.r
    m = c.m
    xᵢ = [_variablesₓᵢ_rm(m, i) for i in 0:m - 1]
    row_matrices = [reduce(_vmult_rm, [xᵢ[i + 1] for i in S], init = ones(Int, 2 ^ m)) for s in 0:r for S in combinations(0:m - 1, s)]
    rows = length(row_matrices)
    cols = length(row_matrices[1])
    G = reshape(vcat(row_matrices...), cols, rows)'
    G = Matrix{Bool}(G)
    return G 
end

function parity_checks(c::ReedMuller)
    H = generator(ReedMuller(c.m - c.r - 1, c.m))
    return H
end

code_n(c::ReedMuller) = 2 ^ c.m
code_k(c::ReedMuller) = sum(binomial.(c.m, 0:c.r))
distance(c::ReedMuller) = 2 ^ (c.m - c.r)
rate(c::ReedMuller) = code_k(c::ReedMuller) / code_n(c::ReedMuller)
