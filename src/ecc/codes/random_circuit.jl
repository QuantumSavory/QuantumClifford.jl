using Random: AbstractRNG, GLOBAL_RNG

struct RandomCircuitCode <: AbstractECC
    arrange::NTuple{N,Int} where {N}
    connect::Union{Val{:alltoall},Val{:brickwork}}
    circ::Vector{QuantumClifford.AbstractOperation}
    encode_qubits::AbstractArray
    # it will be nicer if we can use CartesianIndex for encode_qubits here,
    # but its conversion to LinearIndex is limited, not supporting non-one step.

    function RandomCircuitCode(rng::AbstractRNG, n::Int, connect::Val{:alltoall}, ngates::Int, k::Int)
        new((n,), Val(:alltoall), random_all_to_all_clifford_circuit(rng, n, ngates), collect(1:k))
    end

    function RandomCircuitCode(n::Int, connect::Val{:alltoall}, ngates::Int, k::Int)
        new((n,), Val(:alltoall), random_all_to_all_clifford_circuit(n, ngates), collect(1:k))
    end

    function RandomCircuitCode(rng::AbstractRNG, n::Int, connect::Val{:alltoall}, ngates::Int, encode_qubits::AbstractArray)
        new((n,), Val(:alltoall), random_all_to_all_clifford_circuit(rng, n, ngates), encode_qubits)
    end

    function RandomCircuitCode(n::Int, connect::Val{:alltoall}, ngates::Int, encode_qubits::AbstractArray)
        new((n,), Val(:alltoall), random_all_to_all_clifford_circuit(n, ngates), encode_qubits)
    end

    function RandomCircuitCode(rng::AbstractRNG, arrange::NTuple{N,Int} where {N}, connect::Val{:brickwork}, nlayers::Int, encode_qubits::AbstractArray)
        new(arrange, Val(:brickwork), random_brickwork_clifford_circuit(rng, arrange, nlayers), encode_qubits)
    end

    function RandomCircuitCode(arrange::NTuple{N,Int} where {N}, connect::Val{:brickwork}, nlayers::Int, encode_qubits::AbstractArray)
        new(arrange, Val(:brickwork), random_brickwork_clifford_circuit(arrange, nlayers), encode_qubits)
    end
end

function random_brickwork_clifford_circuit(rng::AbstractRNG, arrange::NTuple{N,Int} where {N}, nlayers::Int)
    circ = QuantumClifford.AbstractOperation[]
    cartesian = CartesianIndices(arrange)
    dim = length(arrange)
    nqubits = prod(arrange)
    for i in 1:nlayers
        gate_direction = (i - 1) % dim + 1
        l = arrange[gate_direction]
        brickwise_parity = dim == 1 ? i % 2 : 1 - (i ÷ dim) % 2
        for j in 1:nqubits
            cardj = collect(cartesian[j].I)
            if cardj[gate_direction] % 2 == brickwise_parity && cardj[gate_direction] != l # open boundary
                cardk = cardj
                cardk[gate_direction] = cardk[gate_direction] + 1
                k = LinearIndices(cartesian)[cardk...]
                push!(circ, SparseGate(random_clifford(rng, 2), [j, k]))
            end
        end
    end
    circ
end

random_brickwork_clifford_circuit(arrange::NTuple{N,Int} where {N}, nlayers::Int) = random_brickwork_clifford_circuit(GLOBAL_RNG, arrange, nlayers)

function random_all_to_all_clifford_circuit(rng::AbstractRNG, nqubits::Int, ngates::Int)
    circ = QuantumClifford.AbstractOperation[]
    for i in 1:ngates
        j = rand(1:nqubits)
        k = rand(1:nqubits-1)
        push!(circ, SparseGate(random_clifford(rng, 2), [j, (j + k - 1) % nqubits + 1]))
    end
    circ
end

random_all_to_all_clifford_circuit(nqubits::Int, ngates::Int) = random_all_to_all_clifford_circuit(GLOBAL_RNG, nqubits, ngates)

function parity_checks(c::RandomCircuitCode)
    n = code_n(c)
    checks = one(Stabilizer, n)[setdiff(1:n, c.encode_qubits)]
    for op in c.circ
        apply!(checks, op)
    end
    checks
end

iscss(::Type{RandomCircuitCode}) = nothing

code_n(c::RandomCircuitCode) = prod(c.arrange)

code_k(c::RandomCircuitCode) = length(c.encode_qubits)
