import Oscar: free_group, cyclic_group, direct_product, small_group_identification, describe, order, gens, quo,
polynomial_ring, matrix, GF, transpose, laurent_polynomial_ring, ideal
function load_oscar_codes()
  #@info "Add group theoretic codes requiring Oscar"
  # [[72, 8, 9]] 2BGA code taken from Table I Block 1 of [lin2024quantum](@cite)
  F = free_group(["r"])
  r = gens(F)[1]
  G, = quo(F, [r^36])
  GA = group_algebra(GF(2), G)
  r = gens(G)[1]
  a = [one(G), r^28]
  b = [one(G), r, r^18, r^12, r^29, r^14]
  t1b1 = twobga_from_fp_group(a, b, GA)

  # [[54, 6, 9]] 2BGA code taken from Table I Block 3 of [lin2024quantum](@cite)
  F = free_group(["r"])
  r = gens(F)[1]
  G, = quo(F, [r^27])
  GA = group_algebra(GF(2), G)
  r = gens(G)[1]
  a = [one(G), r, r^3, r^7]
  b = [one(G), r, r^12, r^19]
  t1b3 = twobga_from_fp_group(a, b, GA)

  # [[16, 4, 4]] 2BGA taken from Appendix C, Table II of [lin2024quantum](@cite)
  F = free_group(["x", "s"])
  x, s = gens(F)
  G, = quo(F, [x^4, s^2, x * s * x^-1 * s^-1])
  GA = group_algebra(GF(2), G)
  x, s = gens(G)
  a = [one(G), x]
  b = [one(G), x, s, x^2, s*x, x^3]
  tb21 = twobga_from_fp_group(a, b, GA)

  # [[32, 8, 4]] 2BGA taken from Appendix C, Table II of [lin2024quantum](@cite)
  F = free_group(["x", "s"])
  x, s = gens(F)
  G, = quo(F, [x^8, s^2, x * s * x^-1 * s^-1])
  GA = group_algebra(GF(2), G)
  x, s = gens(G)
  a = [one(G), x^6]
  b = [one(G), s * x^7, s * x^4, x^6, s * x^5, s * x^2]
  tb22 = twobga_from_fp_group(a, b, GA)

  # Examples of Abelian 2BGA codes constructed from the Direct Product of two cyclic groups, denoted as `C₂ₘ = Cₘ × C₂`.
  # [[56, 8, 7]] 2BGA taken from Appendix C, Table II of [lin2024quantum](@cite)
  m = 14; n = 2
  C₁₄ = cyclic_group(m)
  C₂ = cyclic_group(n)
  G = direct_product(C₁₄, C₂)
  GA = group_algebra(GF(2), G)
  x, s = gens(GA)[1], gens(GA)[3]
  a = [one(GA), x^8]
  b = [one(GA), x^7, s, x^8, x^9, s * x^4]
  dprod1 = twobga_from_direct_product(a, b, GA)

  # [[48, 24, 2]] 2BGA taken from Appendix C, Table II of [lin2024quantum](@cite)
  m = 12; n = 2
  C₁₂ = cyclic_group(m)
  C₂ = cyclic_group(n)
  G = direct_product(C₁₂, C₂)
  GA = group_algebra(GF(2), G)
  x, s = gens(GA)[1], gens(GA)[4]
  a = [one(GA), s * x^6]
  b = [one(GA), x^3, s * x^6, x^4, s * x^9, s * x^10]
  dprod2 = twobga_from_direct_product(a, b, GA)

  # 2BGA codes using non-abelian groups: Table III of [lin2024quantum](@cite)
  # [[24, 8, 3]]
  m = 6
  F = free_group(["r", "s"])
  r, s = gens(F)
  G, = quo(F, [r^m, s^2, (r*s)^2])
  GA = group_algebra(GF(2), G)
  r, s = gens(G)
  a = [one(G), r^4]
  b = [one(G), s*r^4, r^3, r^4, s*r^2, r]
  nonabel1 = twobga_from_fp_group(a, b, GA)

  # [[24, 12, 2]]
  F = free_group(["r", "s"])
  r, s = gens(F)
  G, = quo(F, [r^m, s^2, (r*s)^2])
  GA = group_algebra(GF(2), G)
  r, s = gens(G)
  a = [one(G), r^3]
  b = [one(G), s*r, r^3, r^4, s*r^4, r]
  nonabel2 = twobga_from_fp_group(a, b, GA)

  # [[32, 8, 4]]
  m = 8
  F = free_group(["r", "s"])
  r, s = gens(F)
  G, = quo(F, [r^m, s^2, (r*s)^2])
  GA = group_algebra(GF(2), G)
  r, s = gens(G)
  a = [one(G), r^2]
  b = [one(G), s*r^5, s*r^4, r^2, s*r^7, s*r^6]
  nonabel3 = twobga_from_fp_group(a, b, GA)

  # [[32, 16, 2]]
  F = free_group(["r", "s"])
  r, s = gens(F)
  G, = quo(F, [r^m, s^2, (r*s)^2])
  GA = group_algebra(GF(2), G)
  r, s = gens(G)
  a = [one(G), r^4]
  b = [one(G), s*r^3, s*r^6, r^4, s*r^7, s*r^2]
  nonabel4 = twobga_from_fp_group(a, b, GA)

  append!(test_twobga_codes, [t1b1, t1b3, tb21, tb22, dprod1, dprod2, nonabel1, nonabel2, nonabel3, nonabel4])

  # Homological Product Codes
  # [[117, 9, 4]] from [xu2024fastparallelizablelogicalcomputation](@cite)
  R, x = polynomial_ring(GF(2), "x")
  l₁ = 3
  H₁ = matrix(R, 2, 3, [x^2 x^2 x^2;
                        x   x^2  0])

  # [[225, 9, 6]] from [xu2024fastparallelizablelogicalcomputation](@cite)
  R, x = polynomial_ring(GF(2), "x")
  l₂ = 3
  H₂ = matrix(R, 3, 4, [x^2 x^2 x^2   0;
                        x^2   0 x^2  x^2;
                        x^2 x^2   x  x^2])

  # 3D Homological product code from [Quintavalle_2021](@cite)
  μ = 2; wc = 3; wr = 4
  c = random_Gallager_ldpc(μ, wc, wr)
  H₃ = matrix(GF(2), c)

  # 3D Homological product code from [Quintavalle_2021](@cite)
  δ₄ = matrix(GF(2), parity_matrix(RepCode(3)))

  # Double Homological product codes
  # [[241, 1, 9]] from Table I of https://arxiv.org/pdf/1805.09271
  δ₁ = [1 1 0;
        0 1 1]

  # [[486, 6, 9]] from Table I of https://arxiv.org/pdf/1805.09271
  δ₂ = [1 1 0;
        0 1 1;
        1 0 1]

  # Trivariate Tricycle Codes from [jacob2025singleshotdecodingfaulttolerantgates](@cite)

  # [[36, 3, 3]] from Table III
  F₂ = GF(2)
  l₁, m₁, p₁ = 3, 2, 2
  R, (x, y, z) = polynomial_ring(F₂, [:x, :y, :z])
  I = ideal(R, [x^l₁ - 1, y^m₁ - 1, z^p₁ - 1])
  S, _ = quo(R, I)
  A₁ = S(1 + x*y*z)
  B₁ = S(1 + x^2*z)
  C₁ = S(1 + x)

  # [[48, 3, 4]] from Table III
  l₂, m₂, p₂ = 4, 2, 2
  I = ideal(R, [x^l₂ - 1, y^m₂ - 1, z^p₂ - 1])
  S, _ = quo(R, I)
  A₂ = S(1 + x)
  B₂ = S(1 + x*z)
  C₂ = S(1 + x*y)

  # [[54, 3, 4]] from Table III
  l₃, m₃, p₃ = 3, 3, 2
  I = ideal(R, [x^l₃ - 1, y^m₃ - 1, z^p₃ - 1])
  S, _ = quo(R, I)
  A₃ = S(1 + y*z)
  B₃ = S(1 + x*z)
  C₃ = S(1 + x*y*z)

  # [[108, 6, 2]] from Table IV
  l₄, m₄, p₄ = 4, 3, 3
  I = ideal(R, [x^l₄ - 1, y^m₄ - 1, z^p₄ - 1])
  S, _ = quo(R, I)
  A₄ = S((1 + x^2)*(1 + x*z))
  B₄ = S(1 + x^2*y^2)
  C₄ = S(1 + x^2*y^2*z^2)

  # Generalized Toric Codes from [liang2025generalizedtoriccodestwisted](@cite)
  # [[12, 4, 2]] from Table I of [liang2025generalizedtoriccodestwisted](@cite)
  R, (x,y) = laurent_polynomial_ring(Native.GF(2), [:x, :y])
  f₁ = 1 + x + x*y
  g₁ = 1 + y + x*y
  α1₁ = (0, 3)
  α2₁ = (2, 1)

  # [[14, 6, 2]] from Table I of [liang2025generalizedtoriccodestwisted](@cite)
  f₂ = 1 + x + y
  g₂ = 1 + y + x
  α1₂ = (0, 7)
  α2₂ = (1, 2)

  # [[96, 4, 12]] from Table I of [liang2025generalizedtoriccodestwisted](@cite)
  f₃ = 1 + x + x^-2*y
  g₃ = 1 + y + x*y^-2
  α1₃ = (0, 12)
  α2₃ = (4, 2)

  # [[98, 6, 12]] from Table I of [liang2025generalizedtoriccodestwisted](@cite)
  f₄ = 1 + x + x^-1*y^2
  g₄ = 1 + y + x^-2*y^-1
  α1₄ = (0,  7)
  α2₄ = (7, 0)

  # [[112, 6, 12]] from Table II of [liang2025generalizedtoriccodestwisted](@cite)
  f₅ = 1 + x + x^-1*y^2
  g₅ = 1 + y + x^-2*y^-1
  α1₅ =(0, 7)
  α2₅ =(8, 2)

  # [[114, 4, 14]] from Table II of [liang2025generalizedtoriccodestwisted](@cite)
  f₆ = 1 + x + x^-3*y
  g₆ = 1 + y + x^-5
  α1₆ = (0,  3)
  α2₆ = (19, 1)

  # Bivariate Bicycle codes using polynomial quotient ring
  # [[72, 12, 6]]
  l_bb₁=6; m_bb₁=6
  R_bb₁, (x, y) = polynomial_ring(Native.GF(2), [:x, :y])
  I_bb₁ = ideal(R_bb₁, [x^l-1, y^m-1])
  S_bb₁, _ = quo(R_bb₁, I_bb₁)
  A_bb₁ = S_bb₁(x^3 + y + y^2)
  B_bb₁ = S_bb₁(y^3 + x + x^2)

  # [[90, 8, 10]]
  l_bb₂=15; m_bb₂=3
  R_bb₂, (x, y) = polynomial_ring(Native.GF(2), [:x, :y])
  I_bb₂ = ideal(R_bb₂, [x^l-1, y^m-1])
  S_bb₂, _ = quo(R_bb₂, I_bb₂)
  A_bb₂ = S_bb₂(x^9 + y   + y^2)
  B_bb₂ = S_bb₂(1   + x^2 + x^7)

  # [[108, 8, 10]]
  l_bb₃=9; m_bb₃=6
  R_bb₃, (x, y) = polynomial_ring(Native.GF(2), [:x, :y])
  I_bb₃ = ideal(R_bb₃, [x^l-1, y^m-1])
  S_bb₃, _ = quo(R_bb₃, I_bb₃)
  A_bb₃ = S_bb₃(x^3 + y + y^2)
  B_bb₃ = S_bb₃(y^3 + x + x^2)

  # [[54, 8, 6]]
  l_bb₄=3; m_bb₄=9
  R_bb₄, (x, y) = polynomial_ring(Native.GF(2), [:x, :y])
  I_bb₄ = ideal(R_bb₄, [x^l-1, y^m-1])
  S_bb₄, _ = quo(R_bb₄, I_bb₄)
  A_bb₄ = S_bb₄(1   + y^2 + y^4)
  B_bb₄ = S_bb₄(y^3 + x   + x^2)

  # [[98, 6, 12]]
  l_bb₅=7; m_bb₅=7
  R_bb₅, (x, y) = polynomial_ring(Native.GF(2), [:x, :y])
  I_bb₅ = ideal(R_bb₅, [x^l-1, y^m-1])
  S_bb₅, _ = quo(R_bb₅, I_bb₅)
  A_bb₅ = S_bb₅(x^3 + y^5 + y^6)
  B_bb₅ = S_bb₅(y^2 + x^3 + x^5)

  # Multivariate Multicycle Codes
  # t = 2; Bivariate Bicycle codes
  # [[72, 12, 6]]
  l_mm₁ =6; m_mm₁ = 6
  R_mm₁, (x, y) = polynomial_ring(Native.GF(2), [:x, :y])
  I_mm₁ = ideal(R_mm₁, [x^l_mm₁-1, y^m_mm₁-1])
  S_mm₁, _ = quo(R_mm₁, I_mm₁)
  A_mm₁ = S_mm₁(x^3 + y + y^2)
  B_mm₁ = S_mm₁(y^3 + x + x^2)

  # [[90, 8, 10]]
  l_mm₂ =15; m_mm₂ = 3
  R_mm₂, (x, y) = polynomial_ring(Native.GF(2), [:x, :y])
  I_mm₂ = ideal(R_mm₂, [x^l_mm₂-1, y^m_mm₂-1])
  S_mm₂, _ = quo(R_mm₂, I_mm₂)
  A_mm₂ = S_mm₂(x^9 + y   + y^2)
  B_mm₂ = S_mm₂(1   + x^2 + x^7)

  # t = 3; Trivariate Tricycle codes
  # [[60, 3, 4]]
  ℓ_mm₃, m_mm₃, p_mm₃ = 5, 2, 2
  F2_mm₃ = GF(2)
  R_mm₃, (x, y, z) = polynomial_ring(F2_mm₃, [:x, :y, :z])
  I_mm₃ = ideal(R_mm₃, [x^ℓ_mm₃ - 1, y^m_mm₃ - 1, z^p_mm₃ - 1])
  S_mm₃, _ = quo(R_mm₃, I_mm₃)
  A_mm₃ = S_mm₃(1 + x*z)
  B_mm₃ = S_mm₃(1 + x*y)
  C_mm₃ = S_mm₃(1 + x*y*z)

  # [[90, 3, 5]]
  ℓ_mm₄, m_mm₄, p_mm₄ = 5, 3, 2
  F2_mm₄ = GF(2)
  R_mm₄, (x, y, z) = polynomial_ring(F2_mm₄, [:x, :y, :z])
  I_mm₄ = ideal(R_mm₄, [x^ℓ_mm₄ - 1, y^m_mm₄ - 1, z^p_mm₄ - 1])
  S_mm₄, _ = quo(R_mm₄, I_mm₄)
  A_mm₄ = S_mm₄(1 + x)
  B_mm₄ = S_mm₄(1 + x*y)
  C_mm₄ = S_mm₄(1 + x^2*y^2*z)

  oscar_code_instance_args = Dict(
      :DDimensionalSurface => [(2, 3), (3, 2), (3, 3), (4, 2)],
      :DDimensionalToric => [(2, 3), (3, 2), (3, 3), (4, 2)],
      :GeneralizedToric => [(f₁, g₁, α1₁, α2₁), (f₂, g₂, α1₂, α2₂), (f₃, g₃, α1₃, α2₃), (f₄, g₄, α1₄, α2₄), (f₅, g₅, α1₅, α2₅), (f₆, g₆, α1₆, α2₆)],
      :HomologicalProduct => [([H₁, transpose(H₁)], l₁), ([H₂, transpose(H₂)], l₂), ([H₃, transpose(H₃)],), ([δ₄, δ₄, δ₄],)],
      :DoubleHomologicalProduct => [(δ₁,), (δ₂,)],
      :TrivariateTricycle => [(l₁, m₁, p₁, A₁, B₁, C₁), (l₂, m₂, p₂, A₂, B₂, C₂), (l₃, m₃, p₃, A₃, B₃, C₃), (l₄, m₄, p₄, A₄, B₄, C₄)],
      :BivariateBicycleViaPoly => [(l_bb₁, m_bb₁, A_bb₁, B_bb₁), (l_bb₂, m_bb₂, A_bb₂, B_bb₂), (l_bb₃, m_bb₃, A_bb₃, B_bb₃), (l_bb₄, m_bb₄, A_bb₄, B_bb₄), (l_bb₅, m_bb₅, A_bb₅, B_bb₅)],
      :MultivariateMulticycle =>[([l_mm₁,m_mm₁], [A_mm₁, B_mm₁]), ([l_mm₂,m_mm₂], [A_mm₂, B_mm₂]), ([ℓ_mm₃, m_mm₃, p_mm₃], [A_mm₃, B_mm₃, C_mm₃]), ([ℓ_mm₄, m_mm₄, p_mm₄], [A_mm₄, B_mm₄, C_mm₄])]
  )
  merge!(code_instance_args, oscar_code_instance_args)
  return oscar_code_instance_args
end
const oscar_code_instance_args = load_oscar_codes()
