using Test
using QuantumClifford
using QuantumClifford.ECC
using InteractiveUtils

# generate instances of all implemented codes to make sure nothing skips being checked

# we does not include smaller random circuit code, because some of them has a bad distance and fails the TableDecoder test
const random_circuit_args = repeat([
        (20, Val(:alltoall), 200, 1), (40, Val(:alltoall), 200, [1, 20]),
        ((20,), Val(:brickwork), 50, [1]), ((20,), Val(:brickwork), 50, 1:2:20),
        ((5, 5), Val(:brickwork), 50, [1]), ((3, 3, 3), Val(:brickwork), 50, [1])
    ], 10) # repeat for more randomness


random_circuit_code_args = [map(f -> getfield(random_circuit_code(c...), f), fieldnames(CircuitCode)) for c in random_circuit_args]

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
