using Test
using QuantumClifford
using QuantumClifford.ECC
using InteractiveUtils

# generate instances of all implemented codes to make sure nothing skips being checked

const code_instance_args = Dict(
    Toric => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Surface => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Gottesman => [3, 4, 5],
    CSS => (c -> (parity_checks_x(c), parity_checks_z(c))).([Shor9(), Steane7(), Toric(4,4)]),
    Concat => [(Perfect5(), Perfect5()), (Perfect5(), Steane7()), (Steane7(), Cleve8()), (Toric(2,2), Shor9())],
    RandomCircuitCode => repeat([
            (10, Val(:alltoall), 30, 8), (10, Val(:alltoall), 30, 1:2:7),
            ((7,), Val(:brickwork), 3, [1, 7]), ((20,), Val(:brickwork), 4, 1:4:20),
            ((4, 6), Val(:brickwork), 4, 1:8)
        ], 5), # repeat for more randomness
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
