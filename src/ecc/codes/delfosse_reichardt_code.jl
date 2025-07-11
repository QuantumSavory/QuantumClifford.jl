"""
    $TYPEDEF

The `[[8rp, (8r − 2)p − 2m, 4]]` Delfosse-Reichardt code is derived from the classical
Reed-Muller code and is used to construct quantum stabilizer code.

Delfosse and Reichardt utilize the `[8, 4, 4]` Reed-Muller code to construct `[[8p, 6(p−1), 4]]`
self-dual CSS quantum codes for `p≥2`, and the `[16, 11, 4]` Reed-Muller code to construct
`[[16p, 14p − 8, 4]]` self-dual CSS quantum codes for `p≥1`. 

To improve the generality of the code construction, we proposed that these codes be generalized
using the Reed-Muller code as the base matrix. Rather than hardcoding the `[8, 4, 4]` or `[16, 11, 4]`
Reed-Muller codes, users should be able to input parameters `r` and `m`, thus enhancing the versatility
of the code. This leads to the following generalized code: `[[8rp, (8r − 2)p − 2m, 4]]` Delfosse-Reichardt
code.

The `[[8p, 6(p − 1), 4]]` and `[[16p, 14p − 8, 4]]` codes were introduced by Delfosse and
Reichardt in the paper *Short Shor-style syndrome sequences* [delfosse2020short](@cite). The
parameter `p` specifies the **number of blocks** in the code construction.

# [[8p, 6(p−1), 4]] code family

An `[[16, 6, 4]]` Delfosse-Reichardt code of from `[[8p, 6(p−1), 4]]` code family
from [delfosse2020short](@cite).

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC: DelfosseReichardt, code_n, code_k, parity_checks; # hide

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
julia> using QuantumClifford; using QuantumClifford.ECC: DelfosseReichardt, code_n, code_k, parity_checks; # hide

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
        new(blocks,r,m)
    end
end

function _generalize_delfosse_reichardt_code(blocks::Int, r::Int, m::Int)
    # base matrix: Reed-Muller parity check matrix
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

code_n(c::DelfosseReichardt) = 8*c.blocks*c.r

code_k(c::DelfosseReichardt) = (8*c.r − 2)*c.blocks − 2*c.m

parity_matrix_x(c::DelfosseReichardt) = _generalize_delfosse_reichardt_code(c.blocks, c.r, c.m)

parity_matrix_z(c::DelfosseReichardt) = _generalize_delfosse_reichardt_code(c.blocks, c.r, c.m)
