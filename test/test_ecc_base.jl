using Test
using QuantumClifford
using QuantumClifford.ECC
using InteractiveUtils

import Nemo: GF
import LinearAlgebra
import Hecke: group_algebra, abelian_group

# generate instances of all implemented codes to make sure nothing skips being checked

# We do not include smaller random circuit code because some of them has a bad distance and fails the TableDecoder test
const random_brickwork_circuit_args = repeat([((20,), 50, [1]), ((20,), 50, 1:2:20), ((5, 5), 50, [1]), ((3, 3, 3), 50, [1])], 10)
const random_all_to_all_circuit_args = repeat([(20, 200, 1), (40, 200, [1, 20])], 10)

random_circuit_code_args = vcat(
    [map(f -> getfield(random_brickwork_circuit_code(c...), f), fieldnames(CircuitCode)) for c in random_brickwork_circuit_args],
    [map(f -> getfield(random_all_to_all_circuit_code(c...), f), fieldnames(CircuitCode)) for c in random_all_to_all_circuit_args]
)

# TODO benchmarking them in wiki
# test codes LP04 and LP118 are from https://arxiv.org/pdf/2111.07029
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

# generalized bicyle codes

test_gb_codes = [
    generalized_bicycle_codes([0, 15, 20, 28, 66], [0, 58, 59, 100, 121], 127),
    generalized_bicycle_codes([0, 1, 14, 16, 22], [0, 3, 13, 20, 42], 63),
]

other_lifted_product_codes = []

# from https://arxiv.org/abs/2202.01702v3
l = 63
GA = group_algebra(GF(2), abelian_group(l))
A = zeros(GA, 7, 7)
x = gens(GA)[]
A[LinearAlgebra.diagind(A)] .= x^27
A[LinearAlgebra.diagind(A, -1)] .= x^54
A[LinearAlgebra.diagind(A, 6)] .= x^54
A[LinearAlgebra.diagind(A, -2)] .= GA(1)
A[LinearAlgebra.diagind(A, 5)] .= GA(1)

B = reshape([1 + x + x^6], (1, 1))

push!(other_lifted_product_codes, LPCode(A, B))

const code_instance_args = Dict(
    Toric => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Surface => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Gottesman => [3, 4, 5],
    CSS => (c -> (parity_checks_x(c), parity_checks_z(c))).([Shor9(), Steane7(), Toric(4, 4)]),
    Concat => [(Perfect5(), Perfect5()), (Perfect5(), Steane7()), (Steane7(), Cleve8()), (Toric(2, 2), Shor9())],
    CircuitCode => random_circuit_code_args,
    LPCode => (c -> (c.A, c.B)).(vcat(LP04, LP118, test_gb_codes, other_lifted_product_codes))
)

function all_testablable_code_instances(;maxn=nothing)
    codeinstances = []
    for t in subtypes(QuantumClifford.ECC.AbstractECC)
        for c in get(code_instance_args, t, [])
            codeinstance = t(c...)
            !isnothing(maxn) && nqubits(codeinstance) > maxn && continue
            push!(codeinstances, codeinstance)
        end
    end
    return codeinstances
end
