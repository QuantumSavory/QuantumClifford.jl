"""
    $TYPEDEF

The Delfosse-Reichardt code is derived from the classical *self-orthogonal* Reed-Muller
codes. For parameters `(r,m)` = `(1,3)` and `(2,4)`, [delfosse2020short](@cite) constructs
families of:

- `[[8p, 6(p-1), 4]]` codes requiring `8` measurement rounds.
- `[[16p, 14p-8, 4]]` codes requiring `10` measurement rounds.

Delfosse and Reichardt ([delfosse2020short](@cite)) utilize the `[8, 4, 4]` Reed-Muller code
to construct `[[8p, 6(p−1), 4]]` self-dual CSS quantum codes for `p≥2`, and the `[16, 11, 4]`
Reed-Muller code to construct `[[16p, 14p − 8, 4]]` self-dual CSS quantum codes for `p ≥ 1`. 
The parameter `p` specifies the **number of blocks** in the code construction. To generalize
the code construction, we extended the approach by using self-orthogonal `Reed-Muller` codes as
base matrices for the `Delfosse-Reichardt` code.

!!! note
    Generalization to higher-order RM codes remains an *open problem*, as highlighted in
    [delfosse2020short](@cite): "Find a minimum-length sequence of parity check measurements
    for distance-`7` fault-tolerant error correction with distance-`8` Reed-Muller codes."

# [[8p, 6(p−1), 4]] code family

An `[[16, 6, 4]]` Delfosse-Reichardt code of from `[[8p, 6(p−1), 4]]` code family
from [delfosse2020short](@cite).

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> p = 2; r = 1; m = 3;

julia> c = parity_checks(DelfosseReichardt(p, r, m))
+ XXXXXXXX________
+ ________XXXXXXXX
+ XXXX____XXXX____
+ XX__XX__XX__XX__
+ X_X_X_X_X_X_X_X_
+ ZZZZZZZZ________
+ ________ZZZZZZZZ
+ ZZZZ____ZZZZ____
+ ZZ__ZZ__ZZ__ZZ__
+ Z_Z_Z_Z_Z_Z_Z_Z_

julia> code_n(c), code_k(c)
(16, 6)
```

# [[16p, 14p − 8, 4]] code family

An `[[32, 20, 4]]` Delfosse-Reichardt code of from `[[16p, 14p − 8, 4]]` code family.

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> p = 2; r = 2; m = 4;

julia> c = parity_checks(DelfosseReichardt(p, r, m))
+ XXXXXXXXXXXXXXXX________________
+ ________________XXXXXXXXXXXXXXXX
+ XXXXXXXX________XXXXXXXX________
+ XXXX____XXXX____XXXX____XXXX____
+ XX__XX__XX__XX__XX__XX__XX__XX__
+ X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_
+ ZZZZZZZZZZZZZZZZ________________
+ ________________ZZZZZZZZZZZZZZZZ
+ ZZZZZZZZ________ZZZZZZZZ________
+ ZZZZ____ZZZZ____ZZZZ____ZZZZ____
+ ZZ__ZZ__ZZ__ZZ__ZZ__ZZ__ZZ__ZZ__
+ Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_

julia> code_n(c), code_k(c)
(32, 20)
```

### Fields
    $TYPEDFIELDS
"""
struct DelfosseReichardt <: AbstractCSSCode
    """The number of blocks in the Delfosse-Reichardt CSS code."""
    blocks::Int
    """The order of the classical Reed-Muller code."""
    r::Int
    """The log-length of the classical Reed-Muller code."""
    m::Int
    function DelfosseReichardt(blocks,r,m)
        blocks < 2 && throw(ArgumentError("The number of blocks must be at least 2 to construct a valid code."))
        if r < 0 || r > m
            throw(ArgumentError("Invalid parameters: r must be non-negative and r ≤ m in order to valid code."))
        end
        if !iszero(mod.(parity_matrix(ReedMuller(r,m))*parity_matrix(ReedMuller(r,m))',2))
            throw(ArgumentError("The `Reed-Muller` parity check matrix must be 'self-orthogonal' to construct a self-dual
            CSS `DelfosseReichardt` code. Use `search_self_orthogonal_rm_codes` to search for good parameters for `Reed-Muller` codes
            that provide `self-orthogonal` seeds."))
        end
        new(blocks,r,m)
    end
end

"""
Search for good parameters of `self-orthogonal` Reed-Muller codes.

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; using QECCore: search_self_orthogonal_rm_codes; # hide

julia> search_self_orthogonal_rm_codes(6)
12-element Vector{Tuple{Int64, Int64}}:
 (0, 1)
 (1, 2)
 (1, 3)
 (2, 3)
 (2, 4)
 (3, 4)
 (2, 5)
 (3, 5)
 (4, 5)
 (3, 6)
 (4, 6)
 (5, 6)
```
"""
function search_self_orthogonal_rm_codes(maxₘ::Int)
    good_params = Tuple{Int, Int}[] 
    for m in 1:maxₘ
        for r in 0:m
            try
                RM = ReedMuller(r, m)
                H = parity_matrix(RM)
                if all(iszero, mod.(H*H', 2))
                    push!(good_params, (r, m))
                end
            catch
                continue 
            end
        end
    end
    return good_params
end

function _generalize_delfosse_reichardt_code(blocks::Int, r::Int, m::Int)
    # base matrix: Reed-Muller parity matrix
    H = parity_matrix(ReedMuller(r,m))
    r, c = size(H)
    new_c = blocks*c
    extended_H = zeros(Bool, r+blocks-1, new_c)
    # Create the first 'blocks' rows with 1s for the appropriate block
    for row in 1:blocks
        @inbounds @simd for block in 0:(blocks-1)
            start_c = block*c+1
            end_c = start_c+c-1
            if block == (row - 1)
                extended_H[row, start_c:end_c] .= 1
            end
        end
    end
    # Copy remaining rows directly to each block
    for row in (blocks+1):r+blocks-1
        @inbounds @simd for block in 0:(blocks-1)
            start_c = block*c+1
            end_c = start_c+c-1
            @views extended_H[row, start_c:end_c] .= H[row-blocks+1, :]
        end
    end
    return extended_H
end

function parity_matrix_xz(c::DelfosseReichardt)
    extended_mat = _generalize_delfosse_reichardt_code(c.blocks, c.r, c.m)
    return extended_mat, extended_mat
end

parity_matrix_x(c::DelfosseReichardt) = parity_matrix_xz(c)[1]

parity_matrix_z(c::DelfosseReichardt) = parity_matrix_xz(c)[2]
