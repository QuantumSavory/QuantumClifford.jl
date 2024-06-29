using Random: AbstractRNG, GLOBAL_RNG


"""
`CircuitCode` is defined by a given encoding circuit `circ`.

- `n`: qubit number
- `circ`: the encoding circuit
- `encode_qubits`: the qubits to be encoded
"""
struct CircuitCode <: AbstractECC
    n::Int
    circ::Vector{QuantumClifford.AbstractOperation}
    encode_qubits::AbstractArray
end

iscss(::Type{CircuitCode}) = nothing

code_n(c::CircuitCode) = prod(c.n)

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
Random Clifford circuit code.

The code is generated by a random Clifford circuit `circ` that encode a subset of qubits `encode_qubits` into logical qubits.
The connectivity of the random circuit can be either all-to-all [brown2013short](@cite) or brickwork in some dimensions [gullans2021quantum](@cite).
Each gate in the random circuit is a random 2-qubit Clifford gate.

For `brickwork` connectivity, the qubits are arranged as a lattice, and `lattice_size` contains side length in each dimension.
For example, a 5×5 lattice will have `lattice_size = (5, 5)`.
"""
function random_circuit_code end

function random_circuit_code(rng::AbstractRNG, n::Int, connect::Val{:alltoall}, ngates::Int, k::Int)
    CircuitCode(n, random_all_to_all_clifford_circuit(rng, n, ngates), collect(1:k))
end

function random_circuit_code(n::Int, connect::Val{:alltoall}, ngates::Int, k::Int)
    CircuitCode(n, random_all_to_all_clifford_circuit(n, ngates), collect(1:k))
end

function random_circuit_code(rng::AbstractRNG, n::Int, connect::Val{:alltoall}, ngates::Int, encode_qubits::AbstractArray)
    CircuitCode(n, random_all_to_all_clifford_circuit(rng, n, ngates), encode_qubits)
end

function random_circuit_code(n::Int, connect::Val{:alltoall}, ngates::Int, encode_qubits::AbstractArray)
    CircuitCode(n, random_all_to_all_clifford_circuit(n, ngates), encode_qubits)
end

# TODO it would be nicer if we can use CartesianIndex for `encode_qubits` in brickworks,
# but its conversion to LinearIndex is limited, not supporting non-one step.
function random_circuit_code(rng::AbstractRNG, lattice_size::NTuple{N,Int} where {N}, connect::Val{:brickwork}, nlayers::Int, encode_qubits::AbstractArray)
    CircuitCode(prod(lattice_size), random_brickwork_clifford_circuit(rng, lattice_size, nlayers), encode_qubits)
end

function random_circuit_code(lattice_size::NTuple{N,Int} where {N}, connect::Val{:brickwork}, nlayers::Int, encode_qubits::AbstractArray)
    CircuitCode(prod(lattice_size), random_brickwork_clifford_circuit(lattice_size, nlayers), encode_qubits)
end
