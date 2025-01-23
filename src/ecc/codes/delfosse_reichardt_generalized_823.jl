"""
The `[[8p, 4p − 2, 3]]` Delfosse-Reichardt Generalized `[[8,2,3]]` code is derived from the
quantum `[[8,2,3]]` code. These codes were introduced by Delfosse and Reichardt in the paper
*Short Shor-style syndrome sequences* [delfosse2020short](@cite). The parameter `p` specifies
the **number of blocks** in the code construction.

An `[[16, 6, 3]]` Delfosse-Reichardt Generalized `[[8,2,3]]` code from [delfosse2020short](@cite).

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> p = 2;

julia> c = parity_checks(DelfosseReichardtGeneralized823(p))
+ ZZZZ____________
+ XXXX____________
+ ____ZZZZ________
+ ____XXXX________
+ ________ZZZZ____
+ ________XXXX____
+ ____________ZZZZ
+ ____________XXXX
+ _XYZ_XYZ_XYZ_XYZ
+ _ZXY_ZXY_ZXY_ZXY

julia> code_n(c), code_k(c)
(16, 6)
```
"""
struct DelfosseReichardtGeneralized823 <: AbstractECC
    blocks::Int
    function DelfosseReichardtGeneralized823(blocks)
        blocks < 1 && throw(ArgumentError("The number of blocks must be at least 1 to construct a valid code."))
        new(blocks)
    end
end

function parity_checks(c::DelfosseReichardtGeneralized823)
    H = parity_checks(SmallestColorCode())
    rows, cols = size(H)
    tab = zero(Stabilizer, rows - 2, cols)
    H_rep₁ = parity_checks(SmallestColorCode())[1:4, :]
    H_rep₂ = parity_checks(SmallestColorCode())[5:6, :]
    rows = [hcat(fill(tab, i - 1)..., H_rep₁, fill(tab, c.blocks - i)...) for i in 1:c.blocks]
    D = vcat(rows...)
    E = hcat(fill(H_rep₂, c.blocks)...)
    extended_H = vcat(D, E)
    return extended_H
end

code_n(c::DelfosseReichardtGeneralized823) = 8*c.blocks

code_k(c::DelfosseReichardtGeneralized823) = 4*c.blocks - 2

distance(c::DelfosseReichardtGeneralized823) = 3
