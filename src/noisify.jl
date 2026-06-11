struct NoNoise <: AbstractNoise end

struct CircuitNoise
    single_qubit::AbstractNoise
    two_qubit::AbstractNoise
    idle_noise::AbstractNoise
    measurement::AbstractNoise
    reset::AbstractNoise
end

function CircuitNoise(;
    single_qubit::AbstractNoise = NoNoise(),
    two_qubit::AbstractNoise    = NoNoise(),
    idle_noise::AbstractNoise   = NoNoise(),
    measurement::AbstractNoise  = NoNoise(),
    reset::AbstractNoise = NoNoise()
)
    return CircuitNoise(single_qubit, two_qubit, idle_noise, measurement, reset)
end

CircuitNoise(noise::AbstractNoise) = CircuitNoise(
    single_qubit = noise,
    two_qubit    = noise,
    idle_noise   = noise,
    measurement  = noise,
    reset = noise
)

apply_idle_noise(circuit::AbstractVector, ::NoNoise, nqubits::Integer) = circuit

function apply_idle_noise(circuit::AbstractVector, idle_noise::AbstractNoise, nqubits::Integer)
    isempty(circuit) && return Any[]
    filled_up_to = fill(1, nqubits)
    layers = Dict{Int, Vector{Any}}()
    active_qubits = Dict{Int, Set{Int}}()

    for op in circuit
        if op isa AbstractNoiseOp || op isa VerifyOp || op isa ClassicalXOR
            step = maximum(filled_up_to)
            push!(get!(layers, step, Any[]), op)
            active_qubits[step] = Set(1:nqubits)
            continue
        end

        qs = collect(affectedqubits(op))

        if isempty(qs)
            step = maximum(filled_up_to)
        else
            step = maximum(filled_up_to[qs])
        end

        push!(get!(layers, step, Any[]), op)

        if !isempty(qs)
            union!(get!(active_qubits, step, Set{Int}()), qs)
            filled_up_to[qs] .= step + 1
        end
    end

    output = Any[]
    max_step = maximum(keys(layers))

    for step in 1:max_step
        ops = get(layers, step, Any[])
        active = get(active_qubits, step, Set{Int}())

        idle_qs = setdiff(1:nqubits, collect(active))

        if !isempty(idle_qs)
            push!(output, NoiseOp(idle_noise, idle_qs))
        end

        append!(output, ops)
    end

    return output
end

# code for noisify, two 'versions' essentially, one which just takes in simple noise as input and the other which takes in a CircuitNoise object

noisify(circuit::AbstractVector, noise::AbstractNoise) = reduce(vcat, noisify.(circuit, (noise,)))
noisify(op, noise) = Any[op]
noisify(op::AbstractNoiseOp, noise::AbstractNoise) = Any[op]
noisify(op::ClassicalXOR, noise::AbstractNoise) = Any[op]
noisify(op::VerifyOp, noise::AbstractNoise) = Any[op]

# regular noise

noisify(op::AbstractSingleQubitOperator, noise::AbstractNoise) = Any[NoiseOp(noise, affectedqubits(op)), op]
noisify(op::AbstractTwoQubitOperator, noise::AbstractNoise) = Any[NoiseOp(noise, affectedqubits(op)), op]
noisify(op::AbstractMeasurement, noise::AbstractNoise) = Any[NoiseOp(noise, affectedqubits(op)), op]
noisify(op::Reset, noise::AbstractNoise) = Any[op, NoiseOp(noise, affectedqubits(op))]
noisify(op::AbstractSingleQubitOperator, ::NoNoise) = Any[op]
noisify(op::AbstractTwoQubitOperator, ::NoNoise) = Any[op]
noisify(op::AbstractMeasurement, ::NoNoise) = Any[op]
noisify(op::Reset, ::NoNoise) = Any[op]
# circuit noise


add_noise_op!(out, ::NoNoise, qubits) = out


function add_noise_op!(out, noise::AbstractNoise, qubits)
    push!(out, NoiseOp(noise, qubits))
    return out
end


function noisify(circuit::AbstractVector, noise_model::CircuitNoise; nqubits=nothing)
    if noise_model.idle_noise isa NoNoise
        return reduce(vcat, noisify.(circuit, (noise_model,)))
    end

    nqubits === nothing &&
        error("nqubits must be provided when idle noise is configured")

    idle_noisy_circuit = apply_idle_noise(circuit, noise_model.idle_noise, nqubits)

    return reduce(vcat, noisify.(idle_noisy_circuit, (noise_model,)))
end

noisify(op, ::CircuitNoise) = Any[op]
noisify(op::AbstractNoiseOp, ::CircuitNoise) = Any[op]
noisify(op::ClassicalXOR, ::CircuitNoise) = Any[op]
noisify(op::VerifyOp, ::CircuitNoise) = Any[op]

function noisify(op::AbstractSingleQubitOperator, noise_model::CircuitNoise)
    out = Any[]
    add_noise_op!(out, noise_model.single_qubit, affectedqubits(op))
    push!(out, op)
    return out
end

function noisify(op::AbstractTwoQubitOperator, noise_model::CircuitNoise)
    out = Any[]
    add_noise_op!(out, noise_model.two_qubit, affectedqubits(op))
    push!(out, op)
    return out
end

function noisify(op::AbstractMeasurement, noise_model::CircuitNoise)
    out = Any[]
    add_noise_op!(out, noise_model.measurement, affectedqubits(op))
    push!(out, op)
    return out
end

function noisify(op::Reset, noise_model::CircuitNoise)
    out = Any[]
    push!(out, op)
    add_noise_op!(out, noise_model.reset, affectedqubits(op))
    return out
end
