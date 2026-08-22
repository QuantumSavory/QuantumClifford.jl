"""
    $TYPEDEF

An arbitrary CSS error correction code defined by its `X` and `Z` checks.

### Fields
    $TYPEDFIELDS
"""
struct CSS <: AbstractCSSCode
    """The parity check matrix of the `X` stabilizers."""
    Hx::Matrix{Bool}
    """The parity check matrix of the `Z` stabilizers."""
    Hz::Matrix{Bool}
    function CSS(Hx, Hz)
        owned_Hx = Matrix{Bool}(Hx)
        owned_Hz = Matrix{Bool}(Hz)
        n = size(owned_Hx, 2)
        if n != size(owned_Hz, 2) error("When constructing a CSS quantum code, the two classical codes are required to have the same block size") end
        #if size(Hx,1)+size(Hz,1) >= n error("When constructing a CSS quantum code, the total number of checks (rows) in the parity checks of the two classical codes have to be lower than the block size (the number of columns).") end
        check_allrowscommute(owned_Hx, owned_Hz) || error("The CSS code just created is invalid -- its rows do not commute. This is either a bug in this library, or an inconsistent parity check matrices were provided to the CSS constructor.")
        new(owned_Hx, owned_Hz)
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

function _css_parity_matrix(Hx::Matrix{Bool}, Hz::Matrix{Bool})
    return [Hx falses(size(Hx)); falses(size(Hz)) Hz]
end

parity_matrix(c::CSS) = _css_parity_matrix(c.Hx, c.Hz)

parity_matrix_xz(c::CSS) = (copy(c.Hx), copy(c.Hz))

parity_matrix_x(c::CSS) = copy(c.Hx)
parity_matrix_z(c::CSS) = copy(c.Hz)

code_n(c::CSS) = size(c.Hx,2)
code_s(c::CSS) = size(c.Hx, 1) + size(c.Hz, 1)


# Parity matrix for general CSS codes
function parity_matrix(c::AbstractCSSCode)
    Hx, Hz = parity_matrix_xz(c)
    css = CSS(Hx, Hz)
    return _css_parity_matrix(css.Hx, css.Hz)
end
