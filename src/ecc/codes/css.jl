"""An arbitrary CSS error correcting code defined by its X and Z checks.

```jldoctest
julia> CSS([0 1 1 0; 1 1 0 0], [1 1 1 1]) |> parity_checks
+ _XX_
+ XX__
+ ZZZZ
```"""
struct CSS <: AbstractECC
    Hx::Matrix{Bool}
    Hz::Matrix{Bool}
    function CSS(Hx, Hz)
        n = size(Hx, 2)
        if n != size(Hz, 2) error("When constructing a CSS quantum code, the two classical codes are required to have the same block size") end
        if size(Hx,1)+size(Hz,1) >= n error("When constructing a CSS quantum code, the total number of checks (rows) in the parity checks of the two classical codes have to be lower than the block size (the number of columns).") end
        new(Hx, Hz)
    end
end

function parity_checks(c::CSS)
    extended_Hx = Matrix{Bool}(vcat(c.Hx, zeros(size(c.Hz))))
    extended_Hz = Matrix{Bool}(vcat(zeros(size(c.Hx)), c.Hz))
    Stabilizer(fill(0x0, size(c.Hx, 1) + size(c.Hz, 1)), extended_Hx, extended_Hz)
end

code_n(c::CSS) = size(c.Hx,2)

code_s(c::CSS) = size(c.Hx, 1) + size(c.Hz, 1)
