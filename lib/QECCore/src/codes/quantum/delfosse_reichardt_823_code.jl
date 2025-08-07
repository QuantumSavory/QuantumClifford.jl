"""
    $TYPEDEF

The `[[8p, 4p − 2, 3]]` Delfosse-Reichardt Generalized `[[8,2,3]]` code is derived from the
quantum `[[8,2,3]]` code. These codes were introduced by Delfosse and Reichardt in the paper
*Short Shor-style syndrome sequences* [delfosse2020short](@cite). The parameter `p` specifies
the **number of blocks** in the code construction.

The `[[8, 2, 3]]` non-CSS code serves as the seed code for constructing Delfosse-Reichardt
generalized `[[8p, 4p − 2, 3]]`codes.

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> p = 1;

julia> c = parity_checks(DelfosseReichardt823(p))
+ ZZZZ____
+ XXXX____
+ ____ZZZZ
+ ____XXXX
+ _XYZ_XYZ
+ _ZXY_ZXY

julia> code_n(c), code_k(c)
(8, 2)
```

An `[[16, 6, 3]]` Delfosse-Reichardt Generalized `[[8,2,3]]` code from [delfosse2020short](@cite).

```jldoctest
julia> using QuantumClifford; using QuantumClifford.ECC; # hide

julia> p = 2;

julia> c = parity_checks(DelfosseReichardt823(p))
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

### Fields
    $TYPEDFIELDS
"""
struct DelfosseReichardt823 <: AbstractQECC
    """The number of blocks in the Delfosse-Reichardt generalized [[8, 2, 3]] code."""
    p::Int
    function DelfosseReichardt823(p)
        p < 1 && throw(ArgumentError(THROW_DELFOSSE_823_MIN_BLOCKS))
        new(p)
    end
end

# The `[[8, 2, 3]]` non-CSS code serves as the seed code for constructing Delfosse-Reichardt generalized `[[8p, 4p − 2, 3]]`codes.
_seed₈₂₃ = Bool[0  0  0  0  0  0  0  0  1  1  1  1  0  0  0  0;
                1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0;
                0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1;
                0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0;
                0  1  1  0  0  1  1  0  0  0  1  1  0  0  1  1;
                0  0  1  1  0  0  1  1  0  1  0  1  0  1  0  1]

function parity_matrix(c::DelfosseReichardt823)
    H = _seed₈₂₃
    n = size(H,2)÷2
    Hx = H[:, 1:n]
    Hz = H[:, n+1:2n]
    tab_x = zeros(Bool, size(Hx,1)-2, size(Hx,2))
    tab_z = zeros(Bool, size(Hz,1)-2, size(Hz,2))
    Hx_rep₁ = Hx[1:4, :]
    Hz_rep₁ = Hz[1:4, :]
    Hx_rep₂ = Hx[5:6, :]
    Hz_rep₂ = Hz[5:6, :]
    rows_x = [hcat(fill(tab_x, i-1)..., Hx_rep₁, fill(tab_x, c.p-i)...) for i in 1:c.p]
    rows_z = [hcat(fill(tab_z, i-1)..., Hz_rep₁, fill(tab_z, c.p-i)...) for i in 1:c.p]
    Dx = vcat(rows_x...)
    Dz = vcat(rows_z...)
    D = [Dx Dz]
    Ex = hcat(fill(Hx_rep₂, c.p)...)
    Ez = hcat(fill(Hz_rep₂, c.p)...)
    E = [Ex Ez]
    extended_H = vcat(D, E)
    return extended_H
end

code_n(c::DelfosseReichardt823) = 8*c.p

code_k(c::DelfosseReichardt823) = 4*c.p - 2

distance(c::DelfosseReichardt823) = 3
