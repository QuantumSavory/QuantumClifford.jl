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

insert_idle_noise(circuit::AbstractVector, ::NoNoise) = circuit

function insert_idle_noise(circuit::AbstractVector, idle_noise::AbstractNoise)
    isempty(circuit) && return Any[]
    nqubits = maximum(q for op in circuit for q in affectedqubits(op))
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



noisify(circuit::AbstractVector, noise::AbstractNoise) = reduce(vcat, noisify.(circuit, (noise,)))
noisify(op, noise::AbstractNoise) = Any[op]
noisify(op::AbstractNoiseOp, noise::AbstractNoise) = Any[op]
noisify(op::ClassicalXOR, noise::AbstractNoise) = Any[op]
noisify(op::VerifyOp, noise::AbstractNoise) = Any[op]

noisify(op::AbstractSingleQubitOperator, noise::AbstractNoise) = Any[NoiseOp(noise, affectedqubits(op)), op]
noisify(op::AbstractTwoQubitOperator, noise::AbstractNoise) = Any[NoiseOp(noise, affectedqubits(op)), op]
noisify(op::AbstractMeasurement, noise::AbstractNoise) = Any[NoiseOp(noise, affectedqubits(op)), op]
noisify(op::Reset, noise::AbstractNoise) = Any[op, NoiseOp(noise, affectedqubits(op))]

noisify(op::AbstractSingleQubitOperator, ::NoNoise) = Any[op]
noisify(op::AbstractTwoQubitOperator, ::NoNoise) = Any[op]
noisify(op::AbstractMeasurement, ::NoNoise) = Any[op]
noisify(op::Reset, ::NoNoise) = Any[op]


function noisify(circuit::AbstractVector, noise_model::CircuitNoise)
    if noise_model.idle_noise isa NoNoise
        return reduce(vcat, noisify.(circuit, (noise_model,)))
    end

    idle_noisy_circuit = insert_idle_noise(circuit, noise_model.idle_noise)

    return reduce(vcat, noisify.(idle_noisy_circuit, (noise_model,)))
end

noisify(op, ::CircuitNoise) = Any[op]
noisify(op::AbstractNoiseOp, ::CircuitNoise) = Any[op]
noisify(op::ClassicalXOR, ::CircuitNoise) = Any[op]
noisify(op::VerifyOp, ::CircuitNoise) = Any[op]

noisify(op::AbstractSingleQubitOperator, noise_model::CircuitNoise) = noisify(op, noise_model.single_qubit)
noisify(op::AbstractTwoQubitOperator, noise_model::CircuitNoise) = noisify(op, noise_model.two_qubit)
noisify(op::AbstractMeasurement, noise_model::CircuitNoise) = noisify(op, noise_model.measurement)
noisify(op::Reset, noise_model::CircuitNoise) = noisify(op, noise_model.reset)
