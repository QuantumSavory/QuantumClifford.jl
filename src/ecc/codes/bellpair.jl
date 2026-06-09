"""
`BellPairCode` is a circuit-defined code that treats the input register as `n` physical Bell pairs
(`2n` qubits total) and encodes a chosen set of Bell-pair slots into logical Bell pairs.

The `encode_pairs` field identifies the Bell-pair slots that remain logical. Each slot contributes
two logical qubits, so `code_k(BellPairCode(...)) == 2 * length(encode_pairs)`.

The `encode_pairs` indices are always stored in sorted order; two codes with the same slots in
different order are considered identical.

This is useful for non-convolutional block constructions that operate on Bell-pair registers rather
than on a single stream of qubits.

See also: [`CircuitCode`](@ref), [`random_all_to_all_bellpair_code`](@ref),
[`random_brickwork_bellpair_code`](@ref)
"""
struct BellPairCode <: AbstractQECC
    n::Int
    circ::Vector{QuantumClifford.AbstractOperation}
    encode_pairs::Vector{Int}
    function BellPairCode(n::Int, circ::AbstractVector{<:QuantumClifford.AbstractOperation}, encode_pairs::AbstractArray)
        n < 1 && throw(ArgumentError("n must be positive"))
        circ = Vector{QuantumClifford.AbstractOperation}(circ)
        encode_pairs = sort(collect(encode_pairs))  # normalise order on construction
        all(1 .<= encode_pairs) && all(encode_pairs .<= n) || throw(ArgumentError("encode_pairs must contain Bell-pair slot indices between 1 and n"))
        length(unique(encode_pairs)) == length(encode_pairs) || throw(ArgumentError("encode_pairs must not contain duplicates"))
        new(n, circ, encode_pairs)
    end
end

# implement field-wise equality so struct comparison works correctly
function Base.:(==)(a::BellPairCode, b::BellPairCode)
    a.n == b.n && a.encode_pairs == b.encode_pairs && a.circ == b.circ
end

"""
    BellPairCode(code::AbstractQECC)

Create a `BellPairCode` that encodes `code_k(code)` logical Bell pairs into `code_n(code)` physical
Bell pairs using the encoding circuit of the given `code`.

This uses [`naive_encoding_circuit`](@ref) to generate the circuit, remapped to act on the first
qubit of each physical Bell pair.
"""
function BellPairCode(code::AbstractQECC)
    n = code_n(code)
    k = code_k(code)
    # The naive_encoding_circuit operates on n physical qubits.
    # In BellPairCode, we remap it to act on the "left" qubits (2j-1) of n Bell pairs.
    circ = naive_encoding_circuit(code)
    new_circ = [_remap_qubits(op, j -> 2*j-1) for op in circ]
    # The logical qubits in naive_encoding_circuit are n-k+1:n.
    # So the logical Bell-pair slots are n-k+1:n.
    BellPairCode(n, new_circ, collect(n-k+1:n))
end

# Helper to remap qubits in an operation
function _remap_qubits(op::QuantumClifford.AbstractOperation, f::Function)
    T = typeof(op)
    if T <: QuantumClifford.AbstractSingleQubitOperator
        return T(f(op.q))
    elseif T <: QuantumClifford.AbstractTwoQubitOperator
        return T(f(op.q1), f(op.q2))
    elseif T <: QuantumClifford.sMRZ || T <: QuantumClifford.sMRX || T <: QuantumClifford.sMRY ||
           T <: QuantumClifford.sMZ || T <: QuantumClifford.sMX || T <: QuantumClifford.sMY
        return T(f(op.qubit), op.bit)
    elseif T <: QuantumClifford.sZ || T <: QuantumClifford.sX || T <: QuantumClifford.sY
        return T(f(op.q))
    elseif T <: QuantumClifford.sSWAP
        return T(f(op.q1), f(op.q2))
    else
        # Fallback for other operations that might be used
        fields = fieldnames(T)
        args = Any[]
        for fn in fields
            val = getfield(op, fn)
            if fn in (:q, :q1, :q2, :qubit)
                push!(args, f(val))
            elseif fn == :indices
                push!(args, f.(val))
            else
                push!(args, val)
            end
        end
        return T(args...)
    end
end

# iscss depends on the circuit; return nothing only for the no-circuit case, otherwise we cannot determine CSS-ness statically.
iscss(::Type{BellPairCode}) = nothing

code_n(c::BellPairCode) = 2 * c.n

code_k(c::BellPairCode) = 2 * length(c.encode_pairs)

function _bellpair_qubits(encode_pairs::AbstractVector{Int})
    # use sort (non-mutating) — collect already owns the vector but sort! on
    # a freshly allocated array is misleading style; sort makes intent clearer.
    sort(collect(Iterators.flatten(((2 * i - 1, 2 * i) for i in encode_pairs))))
end

function parity_checks(c::BellPairCode)
    n = c.n
    n2 = 2 * n
    # The initial state is n Bell pairs.
    # Each pair i has stabilizers X_{2i-1}X_{2i} and Z_{2i-1}Z_{2i}.
    full_stabs = zero(Stabilizer, 2*n, 2*n)
    for i in 1:n
        # X_{2i-1}X_{2i} (Row 2i-1)
        full_stabs[2*i-1, 2*i-1] = (true, false)
        full_stabs[2*i-1, 2*i]   = (true, false)
        # Z_{2i-1}Z_{2i} (Row 2i)
        full_stabs[2*i, 2*i-1]   = (false, true)
        full_stabs[2*i, 2*i]     = (false, true)
    end

    # Identify which slots are NOT logical. These provide the parity checks.
    ancilla_slots = setdiff(1:n, c.encode_pairs)
    ancilla_row_indices = Int[]
    for i in ancilla_slots
        push!(ancilla_row_indices, 2*i-1)
        push!(ancilla_row_indices, 2*i)
    end

    # The parity checks are the transformed stabilizers of the ancilla pairs.
    checks = full_stabs[ancilla_row_indices]

    for op in c.circ
        apply!(checks, op)
    end
    checks
end

"""
Random all-to-all Clifford Bell-pair code.

The code of `n` Bell pairs is generated by an all-to-all random Clifford circuit of `ngates` gates
that encodes a subset of Bell-pair slots into logical Bell pairs.

Because of the random picking, the size of `encode_pairs` is the only thing that matters for the
code, referred to as `k` Bell pairs.

See also: [`random_all_to_all_clifford_circuit`](@ref), [`BellPairCode`](@ref)
"""
function random_all_to_all_bellpair_code end

function random_all_to_all_bellpair_code(rng::AbstractRNG, n::Int, ngates::Int, k::Int)
    BellPairCode(n, random_all_to_all_clifford_circuit(rng, 2 * n, ngates), collect(1:k))
end

function random_all_to_all_bellpair_code(n::Int, ngates::Int, k::Int)
    BellPairCode(n, random_all_to_all_clifford_circuit(2 * n, ngates), collect(1:k))
end

function random_all_to_all_bellpair_code(rng::AbstractRNG, n::Int, ngates::Int, encode_pairs::AbstractArray)
    BellPairCode(n, random_all_to_all_clifford_circuit(rng, 2 * n, ngates), encode_pairs)
end

function random_all_to_all_bellpair_code(n::Int, ngates::Int, encode_pairs::AbstractArray)
    BellPairCode(n, random_all_to_all_clifford_circuit(2 * n, ngates), encode_pairs)
end

"""
Random brickwork Clifford Bell-pair code.

The code is generated by a brickwork random Clifford circuit of `nlayers` layers that encodes a
subset of Bell-pair slots into logical Bell pairs.

See also: [`random_brickwork_clifford_circuit`](@ref), [`BellPairCode`](@ref)
"""
function random_brickwork_bellpair_code end

# add k::Int overload to match the all-to-all API shape
function random_brickwork_bellpair_code(rng::AbstractRNG, n::Int, nlayers::Int, k::Int)
    BellPairCode(n, random_brickwork_clifford_circuit(rng, (2 * n,), nlayers), collect(1:k))
end

function random_brickwork_bellpair_code(n::Int, nlayers::Int, k::Int)
    BellPairCode(n, random_brickwork_clifford_circuit((2 * n,), nlayers), collect(1:k))
end

function random_brickwork_bellpair_code(rng::AbstractRNG, n::Int, nlayers::Int, encode_pairs::AbstractArray)
    BellPairCode(n, random_brickwork_clifford_circuit(rng, (2 * n,), nlayers), encode_pairs)
end

function random_brickwork_bellpair_code(n::Int, nlayers::Int, encode_pairs::AbstractArray)
    BellPairCode(n, random_brickwork_clifford_circuit((2 * n,), nlayers), encode_pairs)
end
