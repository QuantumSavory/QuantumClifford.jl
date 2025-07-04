"""
    CSS <: AbstractCSSCode
    CSS(Hx, Hz)

An arbitrary CSS error correction code defined by its X and Z checks.

### Fields
- `Hx::Matrix{Bool}`: The parity check matrix of the X stabilizers.
- `Hz::Matrix{Bool}`: The parity check matrix of the Z stabilizers.
"""
struct CSS <: AbstractCSSCode
    Hx::Matrix{Bool}
    Hz::Matrix{Bool}
    function CSS(Hx, Hz)
        n = size(Hx, 2)
        if n != size(Hz, 2) error("When constructing a CSS quantum code, the two classical codes are required to have the same block size") end
        #if size(Hx,1)+size(Hz,1) >= n error("When constructing a CSS quantum code, the total number of checks (rows) in the parity checks of the two classical codes have to be lower than the block size (the number of columns).") end
        check_allrowscommute(Hx, Hz) || error("The CSS code just created is invalid -- its rows do not commute. This is either a bug in this library, or an inconsistent parity check matrices were provided to the CSS constructor.")
        new(Hx, Hz)
    end
end

function check_allrowscommute(Hx, Hz)
    for rowx in eachrow(Hx)
        for rowz in eachrow(Hz)
            comm = sum(rowx .& rowz)
            isodd(comm) && return false
        end
    end
    return true
end

function parity_matrix(c::CSS)
    extended_Hx = Matrix{Bool}(vcat(c.Hx, zeros(size(c.Hz))))
    extended_Hz = Matrix{Bool}(vcat(zeros(size(c.Hx)), c.Hz))
    return hcat(extended_Hx, extended_Hz)
end

parity_matrix_x(c::CSS) = c.Hx
parity_matrix_z(c::CSS) = c.Hz

code_n(c::CSS) = size(c.Hx,2)
code_s(c::CSS) = size(c.Hx, 1) + size(c.Hz, 1)


# Parity matrix for general CSS codes
parity_matrix(c::AbstractCSSCode) = parity_matrix(CSS(parity_matrix_x(c), parity_matrix_z(c)))
