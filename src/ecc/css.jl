"""Struct for arbitrary CSS error correcting codes.

This struct holds:
    - tab: Boolean matrix with the X part taking up the left side and the Z part taking up the right side"""
struct CSS <: ECC
    tab::Matrix{Bool}
end

function CSS end

"""Creates a CSS code using the two provided matrices where H contains the X checks and G contains the Z checks."""
function CSS(H::Matrix{Bool}, G::Matrix{Bool})::CSS
    Hy, Hx = size(H)
    Gy, Gx = size(G)
    comp_matrix = fill(false, (Hy + Gy, Hx + Gx))
    # comp_matrix = Matrix{Bool}(undef, Hy + Gy, Hx + Gx)
    comp_matrix[1:Hy, 1:Hx] = H
    comp_matrix[Hy+1:end, Hx+1:end] = G
    pcm_stab = Stabilizer(fill(0x0, Hy+Gy), GetXTableau(comp_matrix), GetZTableau(comp_matrix))
    return CSS(comp_matrix)
    # return comp_matrix
end

"""Returns the matrix form of the X and Z checks."""
get_xzs(c::CSS) = c.tab

"""Returns the stabilizer making up the parity check tableau."""
parity_checks(c::CSS) = Stabilizer(fill(0x0, Hy+Gy), GetXTableau(c.tab), GetZTableau(c.tab))

"""Returns the block length of the code."""
code_n(c::CSS) = code_n(Stabilizer(fill(0x0, Hy+Gy), GetXTableau(c.tab), GetZTableau(c.tab)))