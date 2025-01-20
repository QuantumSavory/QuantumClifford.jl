"""
The `[[4p, 2(p − 2), 4]]` Delfosse repetition code is derived from the classical [4, 1, 4]
repetition code and is used to construct quantum stabilizer code.  The [4, 1, 4] repetition
code is a classical error-correcting code where each bit is repeated four times to improve
error detection and correction. It is defined by the following parity-check matrix:

\$\$
\\begin{pmatrix}
1 & 1 & 1 & 1 \\\\
0 & 0 & 1 & 1 \\\\
0 & 1 & 0 & 1
\\end{pmatrix}
\$\$

The `[[4p, 2(p − 2), 4]]` codes were introduced by Delfosse and Reichardt in the paper
*Short Shor-style syndrome sequences* [delfosse2020short](@cite). The parameter `p` represents
the **number of blocks** in the code construction, and it must be a multiple of 2. This is
because the code is constructed by adding *eight qubits at a time*, with each addition
consisting of *two blocks of four qubits*. Since each construction step adds two blocks,
the total number of blocks `p` must always be a multiple of 2 for the code to be valid.
"""
struct DelfosseRepCode <: AbstractECC
    blocks::Int
    function DelfosseRepCode(blocks)
        blocks < 2 && throw(ArgumentError("The number of blocks must be at least 2 to construct a valid code."))
        blocks % 2 != 0 && throw(ArgumentError("The number of blocks must be a multiple of 2."))
        new(blocks)
    end
end

function iscss(::Type{DelfosseRepCode})
    return true
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

function parity_checks(c::DelfosseRepCode)
    extended_mat = _extend_414_repetition_code(c.blocks)
    hx, hz = extended_mat, extended_mat
    code = CSS(hx, hz)
    Stabilizer(code)
end

parity_checks_x(c::DelfosseRepCode) = _extend_414_repetition_code(c.blocks)

parity_checks_z(c::DelfosseRepCode) = _extend_414_repetition_code(c.blocks)
