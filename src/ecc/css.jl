"""An arbitrary CSS error correcting code defined by its X and Z checks."""
struct CSS <: AbstractECC
    tab
end

"""Creates a CSS code using the two provided matrices where H contains the X checks and G contains the Z checks."""
function CSS(H, G)
    Hy, Hx = size(H)
    Gy, Gx = size(G)
    comp_matrix = falses(Hy + Gy, Hx + Gx)
    comp_matrix[1:Hy, 1:Hx] = H
    comp_matrix[Hy+1:end, Hx+1:end] = G
    return CSS(comp_matrix)
end

"""Returns the stabilizer making up the parity check tableau."""
function parity_checks(c::CSS)
    Stabilizer(fill(0x0, size(c.tab, 2)), c.tab[:,1:end÷2], c.tab[:,end÷2+1:end])
end

"""Returns the block length of the code."""
code_n(c::CSS) = size(c.tab,2)÷2
