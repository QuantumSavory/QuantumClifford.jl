using Test
using QuantumClifford
using QuantumClifford.ECC
using InteractiveUtils

# generate instances of all implemented codes to make sure nothing skips being checked

# We do not include smaller random circuit code because some of them has a bad distance and fails the TableDecoder test
const random_brickwork_circuit_args = repeat([((20,), 50, [1]), ((20,), 50, 1:2:20), ((5, 5), 50, [1]), ((3, 3, 3), 50, [1])], 10)
const random_all_to_all_circuit_args = repeat([(20, 200, 1), (40, 200, [1, 20])], 10)

random_circuit_code_args = vcat(
    [map(f -> getfield(random_brickwork_circuit_code(c...), f), fieldnames(CircuitCode)) for c in random_brickwork_circuit_args],
    [map(f -> getfield(random_all_to_all_circuit_code(c...), f), fieldnames(CircuitCode)) for c in random_all_to_all_circuit_args]
)

const code_instance_args = Dict(
    Toric => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Surface => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Gottesman => [3, 4, 5],
    CSS => (c -> (parity_checks_x(c), parity_checks_z(c))).([Shor9(), Steane7(), Toric(4, 4)]),
    Concat => [(Perfect5(), Perfect5()), (Perfect5(), Steane7()), (Steane7(), Cleve8()), (Toric(2, 2), Shor9())],
    CircuitCode => random_circuit_code_args
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
