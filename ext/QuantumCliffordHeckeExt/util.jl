""""
Verifies that the specified group elements constitute a regular F-linear representation
of the group G by checking the properties L(a)L(b) = L(ab) and R(a)R(b) = R(ba) for any
elements a, b in the group algebra F[G].
"""
function check_repr_regular_linear(GA::GroupAlgebra)
    a, b = rand(GA), rand(GA)
    # L(a)L(b) = L(ab)
    L_a = representation_matrix(a)
    L_b = representation_matrix(b)
    L_ab = representation_matrix(a*b)
    # R(a)R(b) = R(ba)
    R_a = representation_matrix(a, :right)
    R_b = representation_matrix(b, :right)
    R_ba = representation_matrix(b*a, :right)
    return L_a * L_b == L_ab && R_a * R_b == R_ba
end
