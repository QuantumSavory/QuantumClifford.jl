using Test
using QuantumClifford
using QuantumClifford.ECC
using QuantumClifford.ECC: check_repr_commutation_relation
using InteractiveUtils

import Nemo: GF
import LinearAlgebra
import Hecke: group_algebra, abelian_group, gens

# generate instances of all implemented codes to make sure nothing skips being checked

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
    generalized_bicycle_codes([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 127), # (A1) [[254, 28, 14≤d≤20]]
    generalized_bicycle_codes([0, 1, 14, 16, 22], [0, 3, 13, 20, 42], 63), # (A2) [[126, 28, 8]]
]

test_hcubic_codes = [
    haah_cubic_codes([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 3)
]

other_lifted_product_codes = []

# [[882, 24, d≤24]] code from (B1) in Appendix B of [panteleev2021degenerate](@cite)
l = 63
GA = group_algebra(GF(2), abelian_group(l))
@test check_repr_commutation_relation(GA) # TODO use this check more pervasively throughout the test suite
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
coprimeBB1 = two_block_group_algebra_codes(A, B)

# [[126,12,10]]
l=7; m=9
GA = group_algebra(GF(2), abelian_group([l*m]))
𝜋 = gens(GA)[1]
A = 1   + 𝜋    + 𝜋^58
B = 𝜋^3 + 𝜋^16 + 𝜋^44
coprimeBB2 = two_block_group_algebra_codes(A, B)

test_coprimeBB_codes = [coprimeBB1, coprimeBB2]

# Multivariate Bicycle codes taken from Table 1 of [voss2024multivariatebicyclecodes](@cite)
# Weight-4 [[144, 2, 12]] MBB code
l=8; m=9
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
z = x*y
A = x^3 + y^7
B = x + y^5
weight4mbb = two_block_group_algebra_codes(A, B)

# Weight-5 [96, 4, 8]] MBB code
l=8; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
z = x*y
A = x^6 + x^3
B = z^5 + x^5 + y
weight5mbb = two_block_group_algebra_codes(A, B)

# Weight-6 [[48, 4, 6]] MBB code
l=4; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
z = x*y
A = x^3 + y^5
B = x + z^5 + y^5 + y^2
weight6mbb = two_block_group_algebra_codes(A, B)

# Weight-7 [[30, 4, 5]] MBB code
l=5; m=3
GA = group_algebra(GF(2), abelian_group([l, m]));
x, y = gens(GA)
z = x*y
A = x^4 + x^2
B = x + x^2 + y + z^2 + z^3
weight7mbb = two_block_group_algebra_codes(A, B)

test_mbb_codes = [weight4mbb, weight5mbb, weight6mbb, weight7mbb]

# Bivariate Bicycle codes
# A [[72, 12, 6]] code from Table 3 of [bravyi2024high](@cite).
l=6; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
A = x^3 + y + y^2
B = y^3 + x + x^2
bb1 = two_block_group_algebra_codes(A,B)

# A [[90, 8, 10]] code from Table 3 of [bravyi2024high](@cite).
l=15; m=3
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
A = x^9 + y   + y^2
B = 1   + x^2 + x^7
bb2 = two_block_group_algebra_codes(A,B)

# A [[360, 12, ≤ 24]]  code from Table 3 of [bravyi2024high](@cite).
l=30; m=6
GA = group_algebra(GF(2), abelian_group([l, m]))
x, y = gens(GA)
A = x^9 + y    + y^2
B = y^3 + x^25 + x^26
bb3 = two_block_group_algebra_codes(A,B)

test_bb_codes = [bb1, bb2, bb3]

const code_instance_args = Dict(
    :Toric => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    :Surface => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    :Gottesman => [3, 4, 5],
    :CSS => (c -> (parity_checks_x(c), parity_checks_z(c))).([Shor9(), Steane7(), Toric(4, 4)]),
    :Concat => [(Perfect5(), Perfect5()), (Perfect5(), Steane7()), (Steane7(), Cleve8()), (Toric(2, 2), Shor9())],
    :CircuitCode => random_circuit_code_args,
    :LPCode => (c -> (c.A, c.B)).(vcat(LP04, LP118, test_gb_codes, test_bb_codes, test_mbb_codes, test_coprimeBB_codes, test_hcubic_codes, other_lifted_product_codes)),
    :QuantumReedMuller => [3, 4, 5],
    :StructuredQLDPC => [(4,3,3), (5,4,4), (6,5,5), (7,6,6)],
    :random_structured_qldpc_code => [(6,4,3), (7,5,3), (8,6,3)]
)

function all_testablable_code_instances(;maxn=nothing)
    codeinstances = []
    i = 1
    for t in subtypes(QuantumClifford.ECC.AbstractECC)
        for c in get(code_instance_args, t.name.name, [])
            codeinstance = t(c...)
            !isnothing(maxn) && nqubits(codeinstance) > maxn && continue
            push!(codeinstances, codeinstance)
            #@show i, t, code_n(codeinstance), code_k(codeinstance), code_s(codeinstance), code_n(codeinstance)-code_k(codeinstance)
            i += 1
        end
    end
    return codeinstances
end
