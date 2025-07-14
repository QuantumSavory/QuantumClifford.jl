"""
    $TYPEDEF

The `[[4p, 2(p − 2), 4]]` Delfosse-Reichardt repetition code is derived from the
classical `[4, 1, 4]` repetition code and is used to construct quantum stabilizer code.

The `[4, 1, 4]` repetition code is a classical error-correcting code where each bit
is repeated four times to improve error detection and correction. It is defined by
the following parity-check matrix:

```math
\\begin{pmatrix}
1 & 1 & 1 & 1 \\\\
0 & 0 & 1 & 1 \\\\
0 & 1 & 0 & 1
\\end{pmatrix}
```

The `[[4p, 2(p − 2), 4]]` codes were introduced by Delfosse and Reichardt in the
paper *Short Shor-style syndrome sequences* [delfosse2020short](@cite). The parameter
`p` specifies the **number of blocks** in the code construction. For the code to be
valid, `p` must be a multiple of 2.

An `[[24, 8, 4]]` Delfosse-Reichardt repetition code from [delfosse2020short](@cite).

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC # hide

julia> p = 6;

julia> c = parity_checks(DelfosseReichardtRepCode(6))
+ XXXX____________________
+ ____XXXX________________
+ ________XXXX____________
+ ____________XXXX________
+ ________________XXXX____
+ ____________________XXXX
+ __XX__XX__XX__XX__XX__XX
+ _X_X_X_X_X_X_X_X_X_X_X_X
+ ZZZZ____________________
+ ____ZZZZ________________
+ ________ZZZZ____________
+ ____________ZZZZ________
+ ________________ZZZZ____
+ ____________________ZZZZ
+ __ZZ__ZZ__ZZ__ZZ__ZZ__ZZ
+ _Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z_Z

julia> code_n(c), code_k(c)
(24, 8)
```

### Fields
    $TYPEDFIELDS
"""
struct DelfosseReichardtRepCode <: AbstractCSSCode
    """The number of blocks in the Delfosse-Reichardt Repetition code."""
    blocks::Int
    function DelfosseReichardtRepCode(blocks)
        blocks < 2 && throw(ArgumentError("The number of blocks must be at least 2 to construct a valid code."))
        blocks % 2 != 0 && throw(ArgumentError("The number of blocks must be a multiple of 2."))
        new(blocks)
    end
end

function _extend_414_repetition_code(blocks::Int)
    n = 4*blocks
    H = zeros(Bool, 2+blocks, n)
    @simd for i in 1:blocks
        H[i, (4*(i-1)+1):(4*i)] .= 1
    end
    @simd for i in 1:blocks
        H[blocks+1, (4*(i-1)+3):(4*(i-1)+4)] .= 1
    end
    @simd for i in 1:blocks
        H[blocks+2, (4*(i-1)+2):(4*(i-1)+4):2] .= 1
    end
    H[end, 1:n] .= (1:n) .% 2 .== 0
    return H
end

function parity_matrix_xz(c::DelfosseReichardtRepCode)
    extended_mat = _extend_414_repetition_code(c.blocks)
    return extended_mat, extended_mat
end

parity_matrix_x(c::DelfosseReichardtRepCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::DelfosseReichardtRepCode) = parity_matrix_xz(c)[2]

code_n(c::DelfosseReichardtRepCode) = 4*c.blocks

code_k(c::DelfosseReichardtRepCode) = 2*(c.blocks - 2)
