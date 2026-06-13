module QuantumCliffordQASMExt

import QuantumClifford:
    AbstractOperation,
    parse_qasm3,
    read_qasm3,
    sCNOT,
    sCPHASE,
    sHadamard,
    sId1,
    sInvPhase,
    sMRZ,
    sMZ,
    sPhase,
    sSWAP,
    sX,
    sY,
    sZ

import Quasar

const _SUPPORTED_QASM_OPERATIONS = "id, x, y, z, h, s, sdg, cx, cz, swap, measure, reset"

mutable struct _BitMapping
    mapping::Dict{String, Vector{Int}}
    count::Int
end
_BitMapping() = _BitMapping(Dict{String, Vector{Int}}(), 0)

function _circuit_instruction(type::String, targets::Vector{Int})
    return (
        type=type,
        arguments=Quasar.InstructionArgument[],
        targets=targets,
        controls=Pair{Int, Int}[],
        exponent=1.0,
    )
end

function _quantumclifford_builtin_gates()
    return Dict{String, Quasar.BuiltinGateDefinition}(
        "id" => Quasar.BuiltinGateDefinition("id", String[], ["q"], _circuit_instruction("id", [0])),
        "x" => Quasar.BuiltinGateDefinition("x", String[], ["q"], _circuit_instruction("x", [0])),
        "y" => Quasar.BuiltinGateDefinition("y", String[], ["q"], _circuit_instruction("y", [0])),
        "z" => Quasar.BuiltinGateDefinition("z", String[], ["q"], _circuit_instruction("z", [0])),
        "h" => Quasar.BuiltinGateDefinition("h", String[], ["q"], _circuit_instruction("h", [0])),
        "s" => Quasar.BuiltinGateDefinition("s", String[], ["q"], _circuit_instruction("s", [0])),
        "sdg" => Quasar.BuiltinGateDefinition("sdg", String[], ["q"], _circuit_instruction("sdg", [0])),
        "cx" => Quasar.BuiltinGateDefinition("cx", String[], ["control", "target"], _circuit_instruction("cx", [0, 1])),
        "cz" => Quasar.BuiltinGateDefinition("cz", String[], ["control", "target"], _circuit_instruction("cz", [0, 1])),
        "swap" => Quasar.BuiltinGateDefinition("swap", String[], ["q1", "q2"], _circuit_instruction("swap", [0, 1])),
    )
end

function _qasm_visitor()
    visitor = Quasar.QasmProgramVisitor()
    empty!(Quasar.gate_defs(visitor))
    merge!(Quasar.gate_defs(visitor), _quantumclifford_builtin_gates())
    return visitor
end

function _maybe_rethrow_qasm_error(err)
    if err isa Quasar.QasmVisitorError
        m = match(r"^gate (.+) not defined!$", err.message)
        if m !== nothing
            gate = m.captures[1]
            throw(Quasar.QasmVisitorError("unsupported or undefined OpenQASM gate `$gate`; QuantumClifford QASM import currently supports $_SUPPORTED_QASM_OPERATIONS."))
        end
    end
    rethrow()
end

function _assignment_parts(expr)
    Quasar.head(expr) == :classical_assignment || return nothing
    op, left_hand_side, right_hand_side = expr.args[1].args
    op == Symbol("=") || return nothing
    Quasar.head(right_hand_side) == :measure || return nothing
    return left_hand_side, right_hand_side
end

function _declare_bits!(bits::_BitMapping, visitor, expr)
    Quasar.head(expr) == :classical_declaration || return
    _, var_type = Quasar.declaration_init(visitor, expr)
    var_type isa Quasar.SizedBitVector || return

    var_name = Quasar.name(expr.args[2])
    bit_size = max(visitor(var_type.size), 1)
    bits.mapping[var_name] = collect(bits.count:(bits.count + bit_size - 1))
    for bit_i in 0:(bit_size - 1)
        bits.mapping["$var_name[$bit_i]"] = [bits.count + bit_i]
    end
    bits.count += bit_size
    return
end

function _bit_indices(raw_indices)
    raw_indices isa AbstractRange && return raw_indices
    raw_indices isa AbstractVector && return raw_indices
    return (raw_indices,)
end

function _evaluate_bits(bits::_BitMapping, visitor, bit_expr)::Vector{Int}
    expr_head = Quasar.head(bit_expr)
    if expr_head == :identifier
        bit_name = Quasar.name(bit_expr)
        haskey(bits.mapping, bit_name) || throw(Quasar.QasmVisitorError("no bit register `$bit_name` is defined."))
        return bits.mapping[bit_name]
    elseif expr_head == :indexed_identifier
        bit_name = Quasar.name(bit_expr)
        haskey(bits.mapping, bit_name) || throw(Quasar.QasmVisitorError("no bit register `$bit_name` is defined."))
        indices = visitor(bit_expr.args[2])
        if indices isa StepRange && indices.step > 0 && indices.stop < indices.start
            indices = StepRange(indices.start, indices.step, length(bits.mapping[bit_name]) - 1)
        end
        out = Int[]
        for raw_index in _bit_indices(indices)
            bit_index = raw_index >= 0 ? raw_index : length(bits.mapping[bit_name]) + raw_index
            key = "$bit_name[$bit_index]"
            haskey(bits.mapping, key) || throw(Quasar.QasmVisitorError("invalid bit index `$raw_index` in `$bit_name`."))
            append!(out, bits.mapping[key])
        end
        return out
    elseif expr_head == :array_literal
        out = Int[]
        for arg in bit_expr.args
            append!(out, _evaluate_bits(bits, visitor, arg))
        end
        return out
    else
        throw(Quasar.QasmVisitorError("unsupported measurement target expression `$expr_head`; expected a bit register or indexed bit register."))
    end
end

function _record_assigned_measurement!(destinations::Vector{Union{Nothing, Int}}, bits::_BitMapping, visitor, expr)
    parts = _assignment_parts(expr)
    parts === nothing && return false
    left_hand_side, right_hand_side = parts
    bit_targets = _evaluate_bits(bits, visitor, left_hand_side)
    qubit_targets = Quasar.evaluate_qubits(visitor, right_hand_side.args[1])
    length(bit_targets) == length(qubit_targets) || throw(Quasar.QasmVisitorError(
        "measurement source and destination sizes must match; got $(length(qubit_targets)) qubit(s) and $(length(bit_targets)) bit(s)."
    ))
    append!(destinations, bit_targets)
    return true
end

function _record_measurement_destinations!(destinations::Vector{Union{Nothing, Int}}, bits::_BitMapping, visitor, expr)
    expr_head = Quasar.head(expr)

    if expr_head == :program || expr_head == :scope
        for child in expr.args
            Quasar.head(child) == :end && return
            _record_measurement_destinations!(destinations, bits, visitor, child)
        end
        return
    elseif expr_head == :classical_declaration
        _declare_bits!(bits, visitor, expr)
        _record_assigned_measurement!(destinations, bits, visitor, expr.args[2])
        visitor(expr)
        return
    elseif expr_head == :classical_assignment
        if _record_assigned_measurement!(destinations, bits, visitor, expr)
            visitor(expr)
            return
        end
    elseif expr_head == :measure
        qubit_targets = Quasar.evaluate_qubits(visitor, expr.args[1])
        append!(destinations, fill(nothing, length(qubit_targets)))
        visitor(expr)
        return
    end

    n_before = length(Quasar.instructions(visitor))
    visitor(expr)
    n_after = length(Quasar.instructions(visitor))
    if n_after > n_before
        for instruction in @view Quasar.instructions(visitor)[(n_before + 1):n_after]
            instruction.type == "measure" && push!(destinations, nothing)
        end
    end
    return
end

function _measurement_destinations(parsed)
    visitor = _qasm_visitor()
    bits = _BitMapping()
    destinations = Union{Nothing, Int}[]
    try
        _record_measurement_destinations!(destinations, bits, visitor, parsed)
    catch err
        _maybe_rethrow_qasm_error(err)
    end
    return destinations
end

function _check_basic_instruction(instruction)
    isempty(instruction.arguments) || throw(Quasar.QasmVisitorError(
        "operation `$(instruction.type)` has parameters, but QuantumClifford QASM import currently supports only the non-parameterized Clifford subset: $_SUPPORTED_QASM_OPERATIONS."
    ))
    isempty(instruction.controls) || throw(Quasar.QasmVisitorError(
        "operation `$(instruction.type)` uses gate controls/modifiers that are not supported by this importer; use explicit `cx` or `cz` gates."
    ))
    instruction.exponent == 1.0 || throw(Quasar.QasmVisitorError(
        "operation `$(instruction.type)` has exponent $(instruction.exponent), but only exponent 1 is supported."
    ))
end

function _target(instruction, n_targets::Int)
    length(instruction.targets) == n_targets || throw(Quasar.QasmVisitorError(
        "operation `$(instruction.type)` expects $n_targets target qubit(s), got $(length(instruction.targets))."
    ))
    return instruction.targets .+ 1
end

function _operation_from_instruction(instruction, measurement_bit)::AbstractOperation
    _check_basic_instruction(instruction)
    type = instruction.type
    if type == "id"
        return sId1(only(_target(instruction, 1)))
    elseif type == "x"
        return sX(only(_target(instruction, 1)))
    elseif type == "y"
        return sY(only(_target(instruction, 1)))
    elseif type == "z"
        return sZ(only(_target(instruction, 1)))
    elseif type == "h"
        return sHadamard(only(_target(instruction, 1)))
    elseif type == "s"
        return sPhase(only(_target(instruction, 1)))
    elseif type == "sdg"
        return sInvPhase(only(_target(instruction, 1)))
    elseif type == "cx"
        control, target = _target(instruction, 2)
        return sCNOT(control, target)
    elseif type == "cz"
        control, target = _target(instruction, 2)
        return sCPHASE(control, target)
    elseif type == "swap"
        q1, q2 = _target(instruction, 2)
        return sSWAP(q1, q2)
    elseif type == "measure"
        qubit = only(_target(instruction, 1))
        return measurement_bit === nothing ? sMZ(qubit) : sMZ(qubit, measurement_bit + 1)
    elseif type == "reset"
        qubit = only(_target(instruction, 1))
        return sMRZ(qubit)
    else
        throw(Quasar.QasmVisitorError("unsupported OpenQASM operation `$(instruction.type)`; QuantumClifford QASM import currently supports $_SUPPORTED_QASM_OPERATIONS."))
    end
end

function _quantumclifford_operations(instructions, measurement_destinations)
    operations = AbstractOperation[]
    measurement_index = 0
    for instruction in instructions
        measurement_bit = nothing
        if instruction.type == "measure"
            measurement_index += 1
            measurement_index <= length(measurement_destinations) || throw(Quasar.QasmVisitorError("internal measurement bookkeeping error: missing measurement destination."))
            measurement_bit = measurement_destinations[measurement_index]
        end
        push!(operations, _operation_from_instruction(instruction, measurement_bit))
    end
    measurement_index == length(measurement_destinations) || throw(Quasar.QasmVisitorError("internal measurement bookkeeping error: unused measurement destination."))
    return operations
end

function parse_qasm3(qasm::AbstractString)
    parsed = Quasar.parse_qasm(String(qasm))
    measurement_destinations = _measurement_destinations(parsed)
    visitor = _qasm_visitor()
    try
        visitor(parsed)
    catch err
        _maybe_rethrow_qasm_error(err)
    end
    return _quantumclifford_operations(Quasar.instructions(visitor), measurement_destinations)
end

read_qasm3(path::AbstractString) = parse_qasm3(read(path, String))

end
