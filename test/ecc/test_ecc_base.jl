using Test
using QuantumClifford.ECC.QECCore
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: check_repr_commutation_relation, check_repr_regular_linear
using InteractiveUtils
using SparseArrays

import Nemo: GF, Native
import LinearAlgebra
import Hecke: group_algebra, abelian_group, gens, quo, one, small_group, polynomial_ring, GF, matrix, quo

# generate instances of all implemented codes to make sure nothing skips being checked

const H1 = sparse(Bool[1 0 1 0; 0 1 0 1; 1 1 0 0]);

const H2 = sparse(Bool[1 1 0;0 1 1]);

# We do not include smaller random circuit code because some of them has a bad distance and fails the TableDecoder test
const random_brickwork_circuit_args = repeat([((20,), 50, [1]), ((20,), 50, 1:2:20), ((5, 5), 50, [1]), ((3, 3, 3), 50, [1])], 10)
const random_all_to_all_circuit_args = repeat([(20, 200, 1), (40, 200, [1, 20])], 10)

random_circuit_code_args = vcat(
    [map(f -> getfield(random_brickwork_circuit_code(c...), f), fieldnames(CircuitCode)) for c in random_brickwork_circuit_args],
    [map(f -> getfield(random_all_to_all_circuit_code(c...), f), fieldnames(CircuitCode)) for c in random_all_to_all_circuit_args]
)

# test codes LP04 and LP118 from Appendix A of [raveendran2022finite](@cite),
B04 = Dict(
    7 => [0 0 0 0; 0 1 2 5; 0 6 3 1],
    9 => [0 0 0 0; 0 1 6 7; 0 4 5 2],
    17 => [0 0 0 0; 0 1 2 11; 0 8 12 13],
    19 => [0 0 0 0; 0 2 6 9; 0 16 7 11]
)

B118 = Dict(
    16 => [0 0 0 0 0; 0 2 4 7 11; 0 3 10 14 15],
    21 => [0 0 0 0 0; 0 4 5 7 17; 0 14 18 12 11],
    30 => [0 0 0 0 0; 0 2 14 24 25; 0 16 11 14 13],
)

LP04 = [LPCode(base_matrix, l .- base_matrix', l) for (l, base_matrix) in B04]
LP118 = [LPCode(base_matrix, l .- base_matrix', l) for (l, base_matrix) in B118]

# generalized bicyle codes from (A1) and (A2) Appendix B of [panteleev2021degenerate](@cite).
test_gb_codes = [
    generalized_bicycle_code_as_2bga([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 127), # (A1) [[254, 28, 14≤d≤20]]
    generalized_bicycle_code_as_2bga([0, 1, 14, 16, 22], [0, 3, 13, 20, 42], 63), # (A2) [[126, 28, 8]]
]

test_hcubic_codes = [
    Haah_cubic_code_as_2bga([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 3),
    Haah_cubic_code_as_2bga(8), # (D) [[1024, 30, 13 ≤ d ≤ 32]] Appendix B of [panteleev2021degenerate](@cite).
]

# honeycomb color codes from [eberhardt2024logical](@cite).
test_honeycomb_color_codes = [
    honeycomb_color_code_as_2bga(6 , 6), honeycomb_color_code_as_2bga(9 , 6),
    honeycomb_color_code_as_2bga(12, 6), honeycomb_color_code_as_2bga(12, 9),
]

# Lifted product codes using non-commutative algebras
# These are just made up, they are not known to
# specifically be good, they are not from a paper.
G = small_group(36,1)
GA = group_algebra(GF(2), G)
r, s  = gens(GA);
A = 1 + r
B = 1 + s + r^6 + s^3*r + s*r^7 + s^3*r^5
nonabel1 = two_block_group_algebra_code(A,B)

G = small_group(48,10)
GA = group_algebra(GF(2), G)
r, s  = gens(GA);
A = 1 + s*r^2
B = 1 + r + s^3 + s^4 + s^2*r^5 + s^4*r^6
nonabel2 = two_block_group_algebra_code(A,B)

G = small_group(40,8)
GA = group_algebra(GF(2), G)
r, s  = gens(GA);
A = 1 + s*r^5 + r^5 + s*r^6
B = 1 + s^2 + r + s^2*r^3
nonabel3 = two_block_group_algebra_code(A,B)

test_nonabelian_codes = [nonabel1, nonabel2, nonabel3]

other_lifted_product_codes = []

# [[882, 24, d≤24]] code from (B1) in Appendix B of [panteleev2021degenerate](@cite)
l = 63
GA = group_algebra(GF(2), abelian_group(l))
@test check_repr_commutation_relation(GA) # TODO use this check more pervasively throughout the test suite
@test check_repr_regular_linear(GA) # TODO use this check more pervasively throughout the test suite
A = zeros(GA, 7, 7)
x = gens(GA)[]
A[LinearAlgebra.diagind(A)] .= x^27
A[LinearAlgebra.diagind(A, -1)] .= x^54
A[LinearAlgebra.diagind(A, 6)] .= x^54
A[LinearAlgebra.diagind(A, -2)] .= GA(1)
A[LinearAlgebra.diagind(A, 5)] .= GA(1)
B = reshape([1 + x + x^6], (1, 1))
push!(other_lifted_product_codes, LPCode(A, B))

# coprime Bivariate Bicycle codes from Table 2 of [wang2024coprime](@cite)
# [[108,12,6]]
l=2; m=27
GA = group_algebra(GF(2), abelian_group([l*m]))
𝜋 = gens(GA)[1]
A = 𝜋^2 + 𝜋^5  + 𝜋^44
B = 𝜋^8 + 𝜋^14 + 𝜋^47
coprimeBB1 = two_block_group_algebra_code(A, B)

# [[126,12,10]]
l=7; m=9
GA = group_algebra(GF(2), abelian_group([l*m]))
𝜋 = gens(GA)[1]
A = 1   + 𝜋    + 𝜋^58
B = 𝜋^3 + 𝜋^16 + 𝜋^44
coprimeBB2 = two_block_group_algebra_code(A, B)

test_coprimeBB_codes = [coprimeBB1, coprimeBB2]

# Multivariate Bicycle codes taken from Table 1 of [voss2024multivariatebicyclecodes](@cite)
# Weight-4 [[144, 2, 12]] MBB code
l=8; m=9
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
z = x*y
A = x^3 + y^7
B = x + y^5
weight4mbb = two_block_group_algebra_code(A, B)

# Weight-5 [96, 4, 8]] MBB code
l=8; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
z = x*y
A = x^6 + x^3
B = z^5 + x^5 + y
weight5mbb = two_block_group_algebra_code(A, B)

# Weight-6 [[48, 4, 6]] MBB code
l=4; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
z = x*y
A = x^3 + y^5
B = x + z^5 + y^5 + y^2
weight6mbb = two_block_group_algebra_code(A, B)

# Weight-7 [[30, 4, 5]] MBB code
l=5; m=3
GA = group_algebra(GF(2), abelian_group([l, m]));
x, y = gens(GA)
z = x*y
A = x^4 + x^2
B = x + x^2 + y + z^2 + z^3
weight7mbb = two_block_group_algebra_code(A, B)

test_mbb_codes = [weight4mbb, weight5mbb, weight6mbb, weight7mbb]

# Bivariate Bicycle codes
# A [[72, 12, 6]] code from Table 3 of [bravyi2024high](@cite).
l=6; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
A = x^3 + y + y^2
B = y^3 + x + x^2
bb1 = two_block_group_algebra_code(A,B)

# A [[90, 8, 10]] code from Table 3 of [bravyi2024high](@cite).
l=15; m=3
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
A = x^9 + y   + y^2
B = 1   + x^2 + x^7
bb2 = two_block_group_algebra_code(A,B)

# A [[360, 12, ≤ 24]]  code from Table 3 of [bravyi2024high](@cite).
l=30; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
A = x^9 + y    + y^2
B = y^3 + x^25 + x^26
bb3 = two_block_group_algebra_code(A,B)

test_bb_codes = [bb1, bb2, bb3]

# Generalized Hypergraph Product codes
# [[882, 24, 18 ≤ d ≤ 24]] from Appendix B of [panteleev2021degenerate](@cite)
F = GF(2)
R, x = polynomial_ring(F, "x")
n_ghp1 = 7
l_ghp1 = 63
S_ghp1, _ =  quo(R, x^l_ghp1 - 1)
A_ghp1 = matrix(S_ghp1, n_ghp1, n_ghp1,
         [x^27  0     0     0     0     1     x^54
          x^54  x^27  0     0     0     0     1
          1     x^54  x^27  0     0     0     0
          0     1     x^54  x^27  0     0     0
          0     0     1     x^54  x^27  0     0
          0     0     0     1     x^54  x^27  0
          0     0     0     0     1     x^54  x^27])
b_ghp1 = S_ghp1(1 + x + x^6)

# [[882, 48, 16]] from Appendix B of [panteleev2021degenerate](@cite)
F = GF(2)
R, x = polynomial_ring(F, "x")
n_ghp2 = 7
l_ghp2 = 63
S_ghp2, _ =  quo(R, x^l_ghp2 - 1)
A_ghp2 = matrix(S_ghp2, n_ghp2, n_ghp2,
         [x^27   0     0     1     x^18  x^27  1
          1      x^27  0     0     1     x^18  x^27
          x^27   1     x^27  0     0     1     x^18
          x^18   x^27  1     x^27  0     0     1
          1      x^18  x^27  1     x^27  0     0
          0      1     x^18  x^27  1     x^27  0
          0      0     1     x^18  x^27  1     x^27])
b_ghp2 = S_ghp2(1 + x + x^6)

# Generalized Bicycle and Extended GB codes from [koukoulekidis2024smallquantumcodesalgebraic](@cite)
R, x = polynomial_ring(GF(2), :x)
l_gb₁ = 6
a_gb₁ = 1 + x^4
b_gb₁ = 1 + x + x^2 + x^4
c_gb₁ = GeneralizedBicycle(a_gb₁, b_gb₁, l_gb₁)
p_gb₁ = one(R)
l_gb₂ = 9
a_gb₂ = 1 + x^2
b_gb₂ = 1 + x^5
c_gb₂ = GeneralizedBicycle(a_gb₂, b_gb₂, l_gb₂)
p_gb₂ = one(R)
l_gb₃ = 10
a_gb₃ = 1 + x
b_gb₃ = 1 + x^6
c_gb₃ = GeneralizedBicycle(a_gb₃, a_gb₃, l_gb₃)
p_gb₃ = one(R) + x

# Add some codes that require Oscar, hence do not work on Windows

test_twobga_codes = []

# La-cross code polynomial
F = GF(2)
R, x = polynomial_ring(F, "x")
h₂ = 1 + x + x^2
h₃ = 1 + x + x^3
h₄ = 1 + x + x^4

# Generalized Bivariate Bicycle Codes
# [[108, 8, 10]] from [bravyi2024high](@cite)
l1 = 9
m1 = 6
A1 = [(:x, 3), (:y, 1), (:y, 2)] # A = x³ + y + y²
B1 = [(:y, 3), (:x, 1), (:x, 2)] # B = y³ + x + x²

# [[90, 8, 10]] from [bravyi2024high](@cite)
l2 = 15
m2 = 3
A2 = [(:x, 9), (:y, 1), (:y, 2)] # A = x⁹ + y + y²
B2 = [(:y, 0), (:x, 2), (:x, 7)] # B = 1 + x² + x⁷

# [[72, 12, 6]] from from [bravyi2024high](@cite)
l3 = 6
m3 = 6
A3 = [(:x, 3), (:y, 1), (:y, 2)] # A = x³ + y + y²
B3 = [(:y, 3), (:x, 1), (:x, 2)] # B = y³ + x + x²

# [[54, 8, 6]] from [wang2024coprime](@cite)
l4 = 3
m4 = 9
A4 = [(:x, 0), (:y, 2), (:y, 4)] # A = 1 + y² + y⁴
B4 = [(:y, 3), (:x, 1), (:x, 2)] # B = y³ + x + x²

# [[98, 6, 12]] from [wang2024coprime](@cite)
l5 = 7
m5 = 7
A5 = [(:x, 3), (:y, 5), (:y, 6)] # A = x³ + y⁵ + y⁶
B5 = [(:y, 2), (:x, 3), (:x, 5)] # A = y² + x³ + x⁵

const code_instance_args = Dict(
    :Toric => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    :Surface => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    :Gottesman => [(3,), (4,), (5,)],
    :CSS => (c -> (parity_matrix_x(c), parity_matrix_z(c))).([Shor9(), Steane7(), Toric(4, 4)]),
    :Concat => [(Perfect5(), Perfect5()), (Perfect5(), Steane7()), (Steane7(), Cleve8()), (Toric(2, 2), Shor9())],
    :CircuitCode => random_circuit_code_args,
    :LPCode => (c -> (c.A, c.B)).(vcat(LP04, LP118, test_gb_codes, test_bb_codes, test_mbb_codes, test_coprimeBB_codes, test_hcubic_codes, test_twobga_codes, test_honeycomb_color_codes, test_nonabelian_codes, other_lifted_product_codes)),
    :QuantumReedMuller => [(3,), (4,), (5,)],
    :Triangular488 => [(3,), (5,), (7,), (9,), (11,)],
    :Triangular666 => [(3,), (5,), (7,), (9,), (11,)],
    :DelfosseReichardt => [(2,1,3), (2,2,4), (4,3,5), (4,3,6)],
    :DelfosseReichardtRep => [(4,), (6,), (8,), (10,)],
    :DelfosseReichardt823 => [(2,), (3,), (4,), (5,)],
    :QuantumTannerGraphProduct => [(H1, H2),(H2, H2), (H1, H1), (H2, H1)],
    :CyclicQuantumTannerGraphProduct => [(2,), (3,), (4,)],
    :LaCross => [(5,h₂,true), (6,h₂,true), (8,h₂,true), (7,h₃,false), (7,h₃,true), (9,h₃,true), (9,h₄,true), (10,h₄,true), (12,h₄,true)],
    :TillichZemor => [(4,3,3), (5,4,4), (6,5,5), (7,6,6)],
    :BivariateBicycleViaCirculantMat => [(l1, m1, A1, B1), (l2, m2, A2,B2), (l3, m3, A3, B3), (l4, m4, A4, B4), (l5, m5, A5, B5)],
    :GeneralizedHyperGraphProduct => [(A_ghp1, b_ghp1, l_ghp1), (A_ghp2, b_ghp2, l_ghp2)],
    :GeneralizedBicycle => [(a_gb₁, b_gb₁, l_gb₁), (a_gb₂, b_gb₂, l_gb₂), (a_gb₃ ,b_gb₃, l_gb₃)],
    :ExtendedGeneralizedBicycle => [(c_gb₁, 2, p_gb₁), (c_gb₂, 3, p_gb₂), (c_gb₃, 4, p_gb₃)]
)

function concretesubtypes(T::DataType)
    concrete = []
    for t in subtypes(T)
        isempty(subtypes(t)) ? push!(concrete, t) : append!(concrete, concretesubtypes(t))
    end
    return concrete
end

function all_testable_code_instances(; maxn=nothing, instance_args=code_instance_args)
    codeinstances = []
    i = 1
    _code_instance_args = copy(instance_args)
    for t in concretesubtypes(QuantumClifford.ECC.AbstractQECC)
        for c in pop!(_code_instance_args, nameof(t), [])
            codeinstance = t(c...)
            !isnothing(maxn) && nqubits(codeinstance) > maxn && continue
            push!(codeinstances, codeinstance)
            #@show i, t, code_n(codeinstance), code_k(codeinstance), code_s(codeinstance), code_n(codeinstance)-code_k(codeinstance)
            i += 1
        end
    end
    @test isempty(_code_instance_args) # if this fails, then some code instances were not tested
    return codeinstances
end
