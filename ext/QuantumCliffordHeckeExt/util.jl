"""
Checks the commutation relation between the left and right representation matrices
for two elements `a` and `b` in the group algebra `ℱ[G]` with a general group `G`.
It verifies the commutation relation that states, `L(a)·R(b) = R(b)·L(a)`. This
property shows that matrices from the left and right representation sets commute
with each other, which is an important property related to the CSS orthogonality
condition.
"""
function check_repr_commutation_relation(GA::GroupAlgebra)
    a, b = rand(GA), rand(GA)
    # Check commutation relation: L(a) R(b) == R(b) L(a)
    L_a = representation_matrix(a)
    R_a = representation_matrix(b, :right)
    return L_a * R_a == R_a * L_a
end
