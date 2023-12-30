"""An arbitrary CSS error correcting code defined by its X and Z checks."""
struct CSS <: AbstractECC
    Hx
    Hz
    """Creates a CSS code using the two provided matrices where Hx contains the X checks and Hz contains the Z checks."""
    function CSS(Hx, Hz)
        n = size(Hx, 2)
        if n != size(Hz, 2) error("When constructing a CSS quantum code, the two classical codes are required to have the same block size") end
        if size(Hx,1)+size(Hz,1) >= n error("When constructing a CSS quantum code, the total number of checks (rows) in the parity checks of the two classical codes have to be lower than the block size (the number of columns).") end
        return new(Hx, Hz)
    end
end

function boolean_tableau(c::CSS)
    Hx_height, Hx_width = size(c.Hx)
    Hz_height, Hz_width = size(x.Hz)
    checks_matrix = falses(Hx_height + Hz_height, Hx_width + Hz_width)
    checks_matrix[1:Hx_height, 1:Hx_width] = c.Hx
    checks_matrix[Hx_height+1:end, Hx_width+1:end] = c.Hz
    return CSS(checks_matrix)
end

"""Returns the stabilizer making up the parity check tableau."""
function parity_checks(c::CSS)
    extended_Hx = Matrix{Bool}(vcat(c.Hx, zeros(size(c.Hz))))
    extended_Hz = Matrix{Bool}(vcat(zeros(size(c.Hx)), c.Hz))
    Stabilizer(fill(0x0, size(c.Hx, 1) + size(c.Hz, 1)), extended_Hx, extended_Hz)
end

"""Returns the block length of the code."""
code_n(c::CSS) = size(c.Hx,2)

"""Returns the depth of the parity check matrix"""
code_m(c::CSS) = size(c.Hx, 1) + size(c.Hz, 1)

"""Returns the number of encoded qubits"""
code_k(c::CSS) = (2 * size(c.Hx,2)) - code_m(c)
