"""A noise model that allows different noise channels to be assigned to single-qubit gates, two-qubit gates, idle qubits, measurements, and reset operations."""
struct CircuitNoise
    single_qubit::Union{AbstractNoise,Nothing}
    two_qubit::Union{AbstractNoise,Nothing}
    idle_noise::Union{AbstractNoise,Nothing}
    measurement::Union{AbstractNoise,Nothing}
    reset::Union{AbstractNoise,Nothing}
end

function CircuitNoise(;
    single_qubit::Union{AbstractNoise,Nothing} = nothing,
    two_qubit::Union{AbstractNoise,Nothing} = nothing,
    idle_noise::Union{AbstractNoise,Nothing} = nothing,
    measurement::Union{AbstractNoise,Nothing} = nothing,
    reset::Union{AbstractNoise,Nothing} = nothing,
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

skip_idling_noise(op) = false
skip_idling_noise(op::VerifyOp) = true
skip_idling_noise(op::ClassicalXOR) = true
skip_idling_noise(op::AbstractNoiseOp) = true

insert_idle_noise(circuit::AbstractVector, ::Nothing) = circuit


function insert_idle_noise(circuit::AbstractVector, idle_noise::AbstractNoise)
    all_qubits = [q for op in circuit for q in affectedqubits(op)]
    isempty(all_qubits) && return copy(circuit)
    nqubits = maximum(all_qubits)
    filled_up_to = fill(1, nqubits)
    op_steps = Int[]
    active_qubits = Dict{Int, Set{Int}}()

    for op in circuit
        qs = collect(affectedqubits(op))

        if skip_idling_noise(op) || isempty(qs)
            push!(op_steps, maximum(filled_up_to))
            continue
        end

        step = maximum(filled_up_to[qs])

        push!(op_steps, step)
        union!(get!(active_qubits, step, Set{Int}()), qs)

        filled_up_to[qs] .= step + 1
    end

    output = []
    emitted = Set{Int}()

    for (op, step) in zip(circuit, op_steps)
        if !skip_idling_noise(op) && !(step in emitted)
            active = get(active_qubits, step, Set{Int}())
            idle = [q for q in 1:nqubits if !(q in active)]

            if !isempty(idle)
                push!(output, NoiseOp(idle_noise, idle))
            end

            push!(emitted, step)
        end

        push!(output, op)
    end

    output
end

"""Construct a noisy version of `circuit` by inserting noise operations according to the specified noise model."""
function noisify(circuit::AbstractVector, noise_model::CircuitNoise)
    idle_noisy_circuit = insert_idle_noise(circuit, noise_model.idle_noise)
    return reduce(vcat, noisify.(idle_noisy_circuit, (noise_model,)))
end

noisify(circuit::AbstractVector, noise::AbstractNoise) = reduce(vcat, noisify.(circuit, (noise,)))

noisify(op, ::Nothing) = [op]
noisify(op, ::AbstractNoise) = [op]

noisify(op::AbstractSingleQubitOperator, noise::AbstractNoise) = [NoiseOp(noise, affectedqubits(op)), op]
noisify(op::AbstractTwoQubitOperator, noise::AbstractNoise) = [NoiseOp(noise, affectedqubits(op)), op]
noisify(op::AbstractMeasurement, noise::AbstractNoise) = [NoiseOp(noise, affectedqubits(op)), op]
noisify(op::Reset, noise::AbstractNoise) = [op, NoiseOp(noise, affectedqubits(op))]


noisify(op, ::CircuitNoise) = [op]
noisify(op::AbstractNoiseOp, ::CircuitNoise) = [op]
noisify(op::ClassicalXOR, ::CircuitNoise) = [op]
noisify(op::VerifyOp, ::CircuitNoise) = [op]

noisify(op::AbstractSingleQubitOperator, noise_model::CircuitNoise) = noisify(op, noise_model.single_qubit)
noisify(op::AbstractTwoQubitOperator, noise_model::CircuitNoise) = noisify(op, noise_model.two_qubit)
noisify(op::AbstractMeasurement, noise_model::CircuitNoise) = noisify(op, noise_model.measurement)
noisify(op::Reset, noise_model::CircuitNoise) = noisify(op, noise_model.reset)
