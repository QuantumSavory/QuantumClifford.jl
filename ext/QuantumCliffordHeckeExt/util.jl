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

"""
The quotient ring ``\\mathbb{F}_2[x]/(x^\\ell - 1)``, consisting of binary polynomials
modulo ``x^\\ell - 1``, is isomorphic to the algebra of ``\\ell \\times \\ell`` circulant
matrices over ``\\mathbb{F}_2``. Thus, a polynomial ``p(x) = \\sum_{i=0}^{\\ell-1} p_i x^i``
can be mapped to a circulant matrix ``P`` whose first row is ``(p_0, p_{\\ell-1}, ..., p_1)``, with
each subsequent row being a cyclic right shift of the previous row [koukoulekidis2024smallquantumcodesalgebraic](@cite).
"""
function circulant_matrix_from_polynomial_ring(l::Int, h::FqPolyRingElem)
    R = parent(h)
    x = gen(R)
    _, proj = residue_ring(R, R(x)^l-1)
    h = proj(h)
    lifted_h = lift(h)
    coeffs = Int[lift(ZZ, coeff(lifted_h, i)) for i in 0:l-1]
    H = zero_matrix(GF(2), l, l)
    for i in 1:l
        for j in 1:l
            H[i, j] = coeffs[mod1(j-i+1, l)]
        end
    end
    H = [Int(lift(ZZ, H[i,j])) for i in 1:nrows(H), j in 1:ncols(H)]
    return H
end
