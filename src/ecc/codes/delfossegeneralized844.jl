"""
The `[[8p, 6(p − 1), 4]]` Delfosse Genneralized `[8,4,4]` code is derived from
the classical `[8, 4, 4]` Reed-Muller code and is used to construct quantum stabilizer
code. The `[[8p, 6(p − 1), 4]]` codes were introduced by Delfosse and Reichardt in the
paper *Short Shor-style syndrome sequences* [delfosse2020short](@cite). The parameter
`p` specifies the **number of blocks** in the code construction.

An `[[16, 6, 4]]` Delfosse Genneralized `[8,4,4]` Reed-Muller code from [delfosse2020short](@cite).

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> c = parity_checks(DelfosseGeneralized844(2))
+ XXXXXXXX________
+ ________XXXXXXXX
+ ____XXXX____XXXX
+ __XX__XX__XX__XX
+ _X_X_X_X_X_X_X_X
+ ZZZZZZZZ________
+ ________ZZZZZZZZ
+ ____ZZZZ____ZZZZ
+ __ZZ__ZZ__ZZ__ZZ
+ _Z_Z_Z_Z_Z_Z_Z_Z

julia> code_n(c), code_k(c)
(16, 6)
```
"""
struct DelfosseGeneralized844 <: AbstractECC
    blocks::Int
    r::Int
    m::Int
    function DelfosseGeneralized844(blocks,r,m)
        blocks < 2 && throw(ArgumentError("The number of blocks must be at least 2 to construct a valid code."))
        if r < 0 || r > m
            throw(ArgumentError("Invalid parameters: r must be non-negative and r ≤ m in order to valid code."))
        end
        new(blocks,r,m)
    end
end

function iscss(::Type{DelfosseGeneralized844})
    return true
end

function _extend_844_code(blocks::Int, r::Int, m::Int)
    # base matrix: Reed-Muller paritycheck matrix
    H = parity_checks(ReedMuller(r,m))
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

function parity_checks(c::DelfosseGeneralized844)
    extended_mat = _extend_844_code(c.blocks, c.r, c.m)
    hx, hz = extended_mat, extended_mat
    code = CSS(hx, hz)
    Stabilizer(code)
end

code_n(c::DelfosseGeneralized844) = 8*c.blocks

code_k(c::DelfosseGeneralized844) = 6*(c.blocks - 1)

parity_checks_x(c::DelfosseGeneralized844) = _extend_844_code(c.blocks)

parity_checks_z(c::DelfosseGeneralized844) = _extend_844_code(c.blocks)
