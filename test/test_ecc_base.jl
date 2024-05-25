using Test
using QuantumClifford
using QuantumClifford.ECC: GottesmanDistance4
using InteractiveUtils

# generate instances of all implemented codes to make sure nothing skips being checked

const code_instance_args = Dict(
    Toric => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Surface => [(3,3), (4,4), (3,6), (4,3), (5,5)],
    Gottesman => [3, 4, 5],
    CSS => (c -> (parity_checks_x(c), parity_checks_z(c))).([Shor9(), Steane7(), Toric(4,4)]),
    GottesmanDistance4 => [4, 5, 6, 7, 8]
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
