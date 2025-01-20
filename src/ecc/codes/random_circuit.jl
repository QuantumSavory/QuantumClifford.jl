using Random: AbstractRNG, GLOBAL_RNG


"""
`CircuitCode` is defined by a given encoding circuit `circ`.

- `n`: qubit number
- `circ`: the encoding circuit
- `encode_qubits`: the qubits to be encoded

See also: [`random_all_to_all_circuit_code`](@ref), [`random_brickwork_circuit_code`](@ref)
"""
struct CircuitCode <: AbstractECC
    n::Int
    circ::Vector{QuantumClifford.AbstractOperation}
    encode_qubits::AbstractArray
end

iscss(::Type{CircuitCode}) = false

code_n(c::CircuitCode) = c.n

code_k(c::CircuitCode) = length(c.encode_qubits)

function parity_checks(c::CircuitCode)
    n = code_n(c)
    checks = one(Stabilizer, n)[setdiff(1:n, c.encode_qubits)]
    for op in c.circ
        apply!(checks, op)
    end
    checks
end

"""
Random all-to-all Clifford circuit code [brown2013short](@cite).

The code of `n` qubits is generated by an all-to-all random Clifford circuit of `ngates` gates that encodes a subset of qubits `encode_qubits` into logical qubits.

Because of the random picking, the size of `encode_qubits` is the only thing that matters for the code, referred to as `k`.

See also: [`random_all_to_all_clifford_circuit`](@ref), [`CircuitCode`](@ref)
"""
function random_all_to_all_circuit_code end

function random_all_to_all_circuit_code(rng::AbstractRNG, n::Int, ngates::Int, k::Int)
    CircuitCode(n, random_all_to_all_clifford_circuit(rng, n, ngates), collect(1:k))
end

function random_all_to_all_circuit_code(n::Int, ngates::Int, k::Int)
    CircuitCode(n, random_all_to_all_clifford_circuit(n, ngates), collect(1:k))
end

function random_all_to_all_circuit_code(rng::AbstractRNG, n::Int, ngates::Int, encode_qubits::AbstractArray)
    CircuitCode(n, random_all_to_all_clifford_circuit(rng, n, ngates), encode_qubits)
end

function random_all_to_all_circuit_code(n::Int, ngates::Int, encode_qubits::AbstractArray)
    CircuitCode(n, random_all_to_all_clifford_circuit(n, ngates), encode_qubits)
end


"""
Random brickwork Clifford circuit code [brown2013short](@cite).

The code is generated by a brickwork random Clifford circuit of `nlayers` layers that encodes a subset of qubits `encode_qubits` into logical qubits.

See also: [`random_brickwork_clifford_circuit`](@ref), [`CircuitCode`](@ref)
"""
function random_brickwork_circuit_code end

# TODO it would be nicer if we can use CartesianIndex for `encode_qubits` in brickworks,
# but its conversion to LinearIndex is limited, not supporting non-one step.
function random_brickwork_circuit_code(rng::AbstractRNG, lattice_size::NTuple{N,Int} where {N}, nlayers::Int, encode_qubits::AbstractArray)
    CircuitCode(prod(lattice_size), random_brickwork_clifford_circuit(rng, lattice_size, nlayers), encode_qubits)
end

function random_brickwork_circuit_code(lattice_size::NTuple{N,Int} where {N}, nlayers::Int, encode_qubits::AbstractArray)
    CircuitCode(prod(lattice_size), random_brickwork_clifford_circuit(lattice_size, nlayers), encode_qubits)
end
