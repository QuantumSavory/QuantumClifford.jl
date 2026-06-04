module QuantumCliffordQuasarExt

import Quasar
using QuantumClifford
import QuantumClifford: parse_qasm3, read_qasm3

# ---------------------------------------------------------------------------
# Quasar AST helpers
# Each node is a QasmExpression with fields:
#   .head :: Symbol   e.g. :program, :gate_call, :qubit_declaration, ...
#   .args :: Vector{Any}   children (QasmExpression or literals)
# ---------------------------------------------------------------------------

head(e::Quasar.QasmExpression) = e.head
args(e::Quasar.QasmExpression) = e.args

# Return the single child of a node (e.g. :identifier -> "q")
only_child(e::Quasar.QasmExpression) = e.args[1]

# Extract the string name from an :identifier node
function identifier_name(e::Quasar.QasmExpression)
    @assert head(e) == :identifier
    return only_child(e)::String
end

# Extract integer value from an :integer_literal node
function integer_value(e::Quasar.QasmExpression)
    @assert head(e) == :integer_literal
    return only_child(e)::Int
end

# Extract (register_name::String, index_0based::Int) from an :indexed_identifier node
function indexed_identifier(e::Quasar.QasmExpression)
    @assert head(e) == :indexed_identifier
    name = identifier_name(e.args[1])
    idx  = integer_value(e.args[2])
    return name, idx
end

# ---------------------------------------------------------------------------
# Gate dispatch tables
# ---------------------------------------------------------------------------

const _SINGLE_QUBIT_GATES = Dict{String,Function}(
    "id"   => q -> sId1(q),
    "i"    => q -> sId1(q),
    "x"    => q -> sX(q),
    "y"    => q -> sY(q),
    "z"    => q -> sZ(q),
    "h"    => q -> sHadamard(q),
    "s"    => q -> sPhase(q),
    "sdg"  => q -> sInvPhase(q),
)

const _TWO_QUBIT_GATES = Dict{String,Function}(
    "cx"   => (a, b) -> sCNOT(a, b),
    "cnot" => (a, b) -> sCNOT(a, b),
    "cz"   => (a, b) -> sCPHASE(a, b),
    "swap" => (a, b) -> sSWAP(a, b),
)

# ---------------------------------------------------------------------------
# Register bookkeeping
# Walk prog.args looking for :qubit_declaration and :classical_declaration
# ---------------------------------------------------------------------------

function _build_register_maps(prog::Quasar.QasmExpression)
    qubit_offset = Dict{String,Int}()
    bit_offset   = Dict{String,Int}()
    qcount = 0
    bcount = 0

    for stmt in prog.args
        if head(stmt) == :qubit_declaration
            # args: [identifier_expr, integer_literal_expr]
            name = identifier_name(stmt.args[1])
            size = integer_value(stmt.args[2])
            qubit_offset[name] = qcount + 1
            qcount += size

        elseif head(stmt) == :classical_declaration
            # args: [classical_type_expr, identifier_expr]
            # The type node wraps a SizedBitVector{N}; we read size from the type
            name = identifier_name(stmt.args[2])
            type_node = stmt.args[1]   # :classical_type
            # type_node.args[1] is a SizedBitVector{N} -- get N via sizeof or type param
            sized_type = type_node.args[1]
            size = _bit_size(sized_type)
            bit_offset[name] = bcount + 1
            bcount += size
        end
    end

    return qubit_offset, bit_offset
end

# Extract the size N from a SizedBitVector{N} or fall back to 1
function _bit_size(t)
    if t isa DataType && t <: Quasar.SizedBitVector
        return t.parameters[1]::Int
    end
    # Fallback: try to read it as a type parameter
    try
        return Int(t.parameters[1])
    catch
        return 1
    end
end

# 0-based QASM index -> 1-based QC index
function _qidx(qubit_offset, name, idx0)
    haskey(qubit_offset, name) || error("Unknown qubit register '$name'.")
    return qubit_offset[name] + idx0
end

function _cidx(bit_offset, name, idx0)
    haskey(bit_offset, name) || error("Unknown classical bit register '$name'.")
    return bit_offset[name] + idx0
end

# ---------------------------------------------------------------------------
# Collect qubit targets from a :qubit_targets node
# Returns Vector of (name, idx0) pairs
# ---------------------------------------------------------------------------

function _qubit_targets(targets_node::Quasar.QasmExpression)
    @assert head(targets_node) == :qubit_targets
    child = targets_node.args[1]
    if head(child) == :indexed_identifier
        return [indexed_identifier(child)]
    elseif head(child) == :array_literal
        return [indexed_identifier(a) for a in child.args]
    else
        error("Unexpected qubit_targets child: $(head(child))")
    end
end

# ---------------------------------------------------------------------------
# Convert a single top-level statement
# ---------------------------------------------------------------------------

function _convert_stmt(stmt::Quasar.QasmExpression, qubit_offset, bit_offset)
    h = head(stmt)

    if h == :gate_call
        # args: [identifier_expr, arguments_expr, qubit_targets_expr]
        gname      = lowercase(identifier_name(stmt.args[1]))
        args_node  = stmt.args[2]   # :arguments (parameters)
        tgts_node  = stmt.args[3]   # :qubit_targets

        # Reject parameterized gates
        if !isempty(args_node.args)
            throw(ArgumentError(
                "Unsupported parameterized gate '$gname'. " *
                "Only parameter-free Clifford gates are supported."
            ))
        end

        qubits = _qubit_targets(tgts_node)

        if length(qubits) == 1
            fn = get(_SINGLE_QUBIT_GATES, gname, nothing)
            fn === nothing && throw(ArgumentError(
                "Unsupported single-qubit gate '$gname'. " *
                "Supported gates: $(sort(collect(keys(_SINGLE_QUBIT_GATES))))."
            ))
            q = _qidx(qubit_offset, qubits[1]...) # (name, idx0)
            return fn(q)

        elseif length(qubits) == 2
            fn = get(_TWO_QUBIT_GATES, gname, nothing)
            fn === nothing && throw(ArgumentError(
                "Unsupported two-qubit gate '$gname'. " *
                "Supported gates: $(sort(collect(keys(_TWO_QUBIT_GATES))))."
            ))
            q1 = _qidx(qubit_offset, qubits[1]...)
            q2 = _qidx(qubit_offset, qubits[2]...)
            return fn(q1, q2)

        else
            throw(ArgumentError(
                "Unsupported $(length(qubits))-qubit gate '$gname'."
            ))
        end

    elseif h == :classical_assignment
        # Measurement: binary_op node with head :(=)
        # args[1]: binary_op
        #   args[1]: :(=)
        #   args[2]: indexed_identifier (classical bit target)
        #   args[3]: measure node
        #     args[1]: indexed_identifier (qubit source)
        binop = stmt.args[1]
        head(binop) == :binary_op || return nothing
        binop.args[1] == :(=)    || return nothing

        c_node   = binop.args[2]   # :indexed_identifier for classical bit
        meas_node = binop.args[3]  # :measure
        head(meas_node) == :measure || return nothing

        q_node = meas_node.args[1]  # :indexed_identifier for qubit

        c_name, c_idx0 = indexed_identifier(c_node)
        q_name, q_idx0 = indexed_identifier(q_node)

        q = _qidx(qubit_offset, q_name, q_idx0)
        c = _cidx(bit_offset,   c_name, c_idx0)
        return sMZ(q, c)

    elseif h == :reset
        # args[1]: indexed_identifier for qubit
        q_name, q_idx0 = indexed_identifier(stmt.args[1])
        q = _qidx(qubit_offset, q_name, q_idx0)
        return sMRZ(q, 0)

    elseif h in (:qubit_declaration, :classical_declaration, :version)
        return nothing   # handled in register pass or ignored

    else
        throw(ArgumentError(
            "Unsupported OpenQASM 3 statement ':$h'. " *
            "Only gate calls, measurements, and reset are currently supported."
        ))
    end
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    parse_qasm3(src::AbstractString) -> Vector{<:AbstractOperation}

Parse an OpenQASM 3 program from string `src` and return the equivalent
QuantumClifford circuit as a `Vector` of symbolic operations.

Requires `Quasar` to be loaded.

# Example
```julia
using QuantumClifford, Quasar
circuit = parse_qasm3(\"\"\"
    OPENQASM 3.0;
    qubit[2] q;
    bit[2] c;
    h q[0];
    cx q[0], q[1];
    c[0] = measure q[0];
    c[1] = measure q[1];
    \"\"\")
# => [sHadamard(1), sCNOT(1, 2), sMZ(1, 1), sMZ(2, 2)]
```
"""
function QuantumClifford.parse_qasm3(src::AbstractString)
    prog = Quasar.parse_qasm(src)
    qubit_offset, bit_offset = _build_register_maps(prog)
    ops = Any[]
    for stmt in prog.args
        op = _convert_stmt(stmt, qubit_offset, bit_offset)
        op === nothing || push!(ops, op)
    end
    return ops
end

"""
    read_qasm3(path::AbstractString) -> Vector{<:AbstractOperation}

Read an OpenQASM 3 file and return the equivalent QuantumClifford circuit.

Requires `Quasar` to be loaded.

# Example
```julia
using QuantumClifford, Quasar
circuit = read_qasm3("bell.qasm")
frames  = pftrajectories(circuit; trajectories=10_000)
```
"""
function QuantumClifford.read_qasm3(path::AbstractString)
    return QuantumClifford.parse_qasm3(read(path, String))
end

end # module