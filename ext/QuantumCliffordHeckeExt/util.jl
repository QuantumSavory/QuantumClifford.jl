"""
Checks the commutation relation between the left and right representation matrices
for two randomly-sampled elements `a` and `b` in the group algebra `ℱ[G]` with a general group `G`.
It verifies the commutation relation that states, `L(a)·R(b) = R(b)·L(a)`. This
property shows that matrices from the left and right representation sets commute
with each other, which is an important property related to the CSS orthogonality
condition. By default, we use `Hecke`’s `representation_matrix`, though a *custom*
`repr` function is also supported.

"""
function check_repr_commutation_relation(GA::GroupAlgebra; repr::Function=representation_matrix)
    a, b = rand(GA), rand(GA)
    # Check commutation relation: L(a)R(b) = R(b)L(a)
    L_a = repr(a)
    R_b = repr(b, :right)
    return L_a * R_b == R_b * L_a
end

"""
Verifies that the specified group elements constitute a regular *F-linear* representation
of the group `G` by checking the properties `L(a)L(b) = L(ab)` and `R(a)R(b) = R(ba)` for any
elements `a`, `b` in the group algebra `F[G]`. By default, we use `Hecke`’s `representation_matrix`,
though a *custom* `repr` function is also supported.
"""
function check_repr_regular_linear(GA::GroupAlgebra; repr::Function=representation_matrix)
    a, b = rand(GA), rand(GA)
    # L(a)L(b) = L(ab)
    L_a = repr(a)
    L_b = repr(b)
    L_ab = repr(a*b)
    # R(a)R(b) = R(ba)
    R_a = repr(a, :right)
    R_b = repr(b, :right)
    R_ba = repr(b*a, :right)
    return L_a*L_b == L_ab && R_a*R_b == R_ba
end
