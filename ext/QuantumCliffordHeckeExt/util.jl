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

"""
    $TYPEDEF

Generate a sparse *quasi-cyclic* polynomial matrix A where each entry is either 0 or a monomial
``x^i (0 \\leq i < \\ell)``. The matrix is constructed such that each row and column
has at most w non-zero entries. This ensures the resulting parity-check matrices ``H_X``
and ``H_Z`` are ``(w + deg b(x))``-limited. Thus, the Tanner graphs 𝒯_X and 𝒯_Z have girth 6
(no 4-cycles), improving belief propagation decoder performance [panteleev2021degenerate](@cite).
"""
function random_qc_ghp_code_matrix_A(S, b, n::Int, w::Int, ℓ::Int; min_k::Int=10, max_attempts=1000, rng::AbstractRNG=default_rng())
    R = base_ring(S)
    x = gen(R)
    b_lifted = lift(b)
    xₗ₋₁ = R(x)^ℓ - 1
    g = gcd(b_lifted, xₗ₋₁)
    k_b = degree(g)
    k_b == 0 && min_k > 0 && error("The chosen polynomial b(x) yields k_b=0. It cannot generate a code with min_k=$min_k. Choose a different b(x).")
    F, _ = residue_ring(R, g)
    for attempt in 1:max_attempts
        A = _matrix_A(S, n, w, ℓ, rng)
        _meets_ghp_constraints(A, w) || continue
        A_proj = matrix(F, n, n, [F(lift((A[i,j]))) for i in 1:n, j in 1:n])
        rk_A_proj = rank(A_proj)
        k_A = n - rk_A_proj
        k = 2*k_A*k_b
        k >= min_k && return A
    end
    error("Failed to generate a valid QC matrix with k >= $min_k after $max_attempts attempts. Try increasing max_attempts, decreasing min_k, or choosing a different b(x) with higher k_b.")
end

function _matrix_A(S, n, w, ℓ, rng)
    R = base_ring(S)
    x = gen(R)
    temp = [S(0) for _ in 1:n]
    positions = randperm(rng, n)[1:w]
    for pos in positions
        exp = rand(rng, 0:ℓ-1)
        temp[pos] = S(x^exp)
    end
    A = zero_matrix(S, n, n)
    for i in 1:n
        for j in 1:n
            shift_idx = mod(j - i, n) + 1
            A[i, j] = temp[shift_idx]
        end
    end
    return A
end

function _meets_ghp_constraints(A, w::Int)
    n = size(A, 1)
    # Each row has exactly w non-zero entries.
    for i in 1:n
        row_cnt = count(j -> !iszero(A[i, j]), 1:n)
        row_cnt ≠ w && return false
    end
    # Each column has exactly w non-zero entries.
    for j in 1:n
        col_cnt = count(i -> !iszero(A[i, j]), 1:n)
        col_cnt ≠ w && return false
    end
    # All non-zero entries are monomials x^i.
    for i in 1:n, j in 1:n
        if !iszero(A[i, j])
            p = lift(A[i, j])
            if count(!iszero, coefficients(p)) ≠ 1
                return false
            end
        end
    end
    return true
end
