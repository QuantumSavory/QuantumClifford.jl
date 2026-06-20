module ECC

using QECCore
using QuantumClifford


"""
Parity check tableau of a code.
"""
function parity_checks end

"""
Generator matrix of a code.
"""
function generator end

"""
Check whether a code is CSS.
"""
iscss(::Type{<:AbstractECC}) = false
iscss(::Type{<:AbstractCSSCode}) = true
iscss(c::AbstractECC) = iscss(typeof(c))



parity_checks(s::Stabilizer)::Stabilizer = s
parity_checks(c::AbstractECC)::Stabilizer = Stabilizer(parity_matrix(c))

Stabilizer(c::AbstractECC) = parity_checks(c)

MixedDestabilizer(c::AbstractECC; kwarg...) =
    MixedDestabilizer(Stabilizer(c); kwarg...)


nqubits(c::AbstractECC) = code_n(c)

code_n(c::AbstractECC) = code_n(parity_checks(c))
code_n(s::Stabilizer) = nqubits(s)

code_s(c::AbstractECC) = code_s(parity_checks(c))
code_s(s::Stabilizer) = length(s)

"""
Number of logical qubits.
"""
function code_k(stab::Stabilizer)

    _, _, rank_stab = canonicalize!(copy(stab), ranks=true)

    return code_n(stab) - rank_stab
end

code_k(c::AbstractECC) = code_k(parity_checks(c))

rate(c::AbstractECC) = code_k(c) / code_n(c)


function logx_ops(code)
    md = MixedDestabilizer(parity_checks(code))
    logicalxview(md)
end

function logz_ops(code)
    md = MixedDestabilizer(parity_checks(code))
    logicalzview(md)
end



function parity_matrix(code::AbstractECC)
    stab_to_gf2(parity_checks(code))
end



"""
Error-to-logical observable map.
"""
function faults_matrix(stab::Stabilizer)

    md = MixedDestabilizer(stab)

    _, n = size(stab)

    rank_stab = rank(md)
    k = n - rank_stab

    O = falses(2k, 2n)

    logicals = vcat(
        logicalxview(md),
        logicalzview(md)
    )

    x_errors = one(Stabilizer, n; basis=:X)
    z_errors = one(Stabilizer, n)

    errors = vcat(x_errors, z_errors)

    for (i, logical) in enumerate(logicals)
        O[i, :] .= comm(logical::PauliOperator, errors)
    end

    return O
end

faults_matrix(code::AbstractECC) =
    faults_matrix(parity_checks(code))



function pauli_error(n, bits, op)

    p = zero(PauliOperator, n)

    for idx in bits
        p[idx] = op
    end

    return p
end

function isdegenerate(
    stab::Stabilizer,
    errors
)

    syndromes = comm.(Ref(stab), errors)

    return length(Set(syndromes)) != length(errors)
end

function isdegenerate(
    stab::Stabilizer,
    d::Int = 1
)

    n = nqubits(stab)

    errors = [
        pauli_error(n, bits, op)
        for bits in combinations(1:n, d)
        for op in ((true,false), (false,true))
    ]

    return isdegenerate(stab, errors)
end

isdegenerate(code::AbstractECC, args...) =
    isdegenerate(parity_checks(code), args...)

isdegenerate(stab::AbstractStabilizer, args...) =
    isdegenerate(stabilizerview(stab), args...)


# Core functionality
include("circuits.jl")
include("decoder_pipeline.jl")

# Utilities
include("codes/util.jl")

# Code constructions
include("codes/concat.jl")
include("codes/random_circuit.jl")
include("codes/classical/bch.jl")

# qLDPC
include("codes/classical/lifted.jl")
include("codes/qeccs_from_extensions.jl")

# Decoders
include("decoder_correction_gate.jl")

end # module ECC
