using QuantumClifford
using QuantumClifford: AbstractOperation
using Quasar: Quasar, parse_qasm, AbstractVisitor, Qubit, QasmExpression, 
    head, name, ClassicalVariable, declaration_init, AbstractGateDefinition, 
    process_gate_targets, qubit_defs, gate_defs, qubit_mapping, qubit_count, 
    instructions, hasgate, splat_gate_targets, evaluate_qubits, 
    QasmVisitorError, SizedBitVector, evaluate_unary_op

struct CliffordGateDefinition <: AbstractGateDefinition
    qubit_targets::Vector{String}
    body::DataType
end

builtin_gates() = Dict(
    "id" => CliffordGateDefinition(["q"], sId1),
    "x" => CliffordGateDefinition(["q"], sX),
    "y" => CliffordGateDefinition(["q"], sY),
    "z" => CliffordGateDefinition(["q"], sZ),
    "h" => CliffordGateDefinition(["q"], sHadamard),
    "s" => CliffordGateDefinition(["q"], sPhase),
    "sdg" => CliffordGateDefinition(["q"], sInvPhase),
    "cx" => CliffordGateDefinition(["c", "t"], sCNOT),
    "cz" => CliffordGateDefinition(["c", "t"], sCPHASE),
    "swap" => CliffordGateDefinition(["q1", "q2"], sSWAP),
)

mutable struct Qasm2CliffordVisitor <: AbstractVisitor
    bit_mapping::Dict{String, Vector{Int}}
    bit_count::Int
    gate_defs::Dict{String, CliffordGateDefinition}
    qubit_defs::Dict{String, Qubit}
    qubit_mapping::Dict{String, Vector{Int}}
    qubit_count::Int
    instructions::Vector{AbstractOperation}
    function Qasm2CliffordVisitor()
        new(
            Dict{String, Vector{Int}}(),
            0,
            builtin_gates(),
            Dict{String, Qubit}(),
            Dict{String, Vector{Int}}(),
            0,
            AbstractOperation[],
           )
    end
end

bit_mapping(v::Qasm2CliffordVisitor)  = v.bit_mapping
bit_count(v::Qasm2CliffordVisitor)  = v.bit_count
Quasar.gate_defs(v::Qasm2CliffordVisitor)  = v.gate_defs
Quasar.hasgate(v::Qasm2CliffordVisitor, gate_name::String) = haskey(gate_defs(v), gate_name)
Quasar.qubit_defs(v::Qasm2CliffordVisitor) = v.qubit_defs
Quasar.qubit_mapping(v::Qasm2CliffordVisitor)  = v.qubit_mapping
Quasar.qubit_count(v::Qasm2CliffordVisitor)  = v.qubit_count
Quasar.instructions(v::Qasm2CliffordVisitor)  = v.instructions
Base.push!(v::Qasm2CliffordVisitor, instr::AbstractOperation) = push!(instructions(v), instr)


function _evaluate_bits(::Val{:identifier}, v, bit_expr::QasmExpression)
    bit_name = name(bit_expr)
    mapping    = bit_mapping(v)
    haskey(mapping, bit_name) || throw(QasmVisitorError("Missing input variable '$bit_name'.", "NameError"))
    return mapping[bit_name]
end
function _evaluate_bits(::Val{:indexed_identifier}, v, bit_expr::QasmExpression)
    bit_name = name(bit_expr)
    mapping    = bit_mapping(v)
    haskey(mapping, bit_name) || throw(QasmVisitorError("Missing input variable '$bit_name'.", "NameError"))
    bit_ix   = v(bit_expr.args[2])
    bits     = Iterators.flatmap(bit_ix) do rq
        if rq >= 0
            haskey(mapping, bit_name * "[$rq]") || throw(QasmVisitorError("Invalid bit index '$rq' in '$bit_name'.", "IndexError"))
            return mapping[bit_name * "[$rq]"]
        else
            bit_size = length(mapping[bit_name])
            haskey(mapping, bit_name * "[$(bit_size + rq)]") || throw(QasmVisitorError("Invalid bit index '$rq' in '$bit_name'.", "IndexError"))
            return mapping[bit_name * "[$(bit_size + rq)]"]
        end
    end
    return collect(bits)
end
_evaluate_bits(val, v, bit_expr) = throw(QasmVisitorError("unable to evaluate bits for expression $bit_expr."))
function evaluate_bits(v::AbstractVisitor, bit_targets::Vector)
    final_bits = Iterators.map(bit_expr->_evaluate_bits(Val(head(bit_expr)), v, bit_expr), bit_targets)
    return vcat(final_bits...)
end
evaluate_bits(v::AbstractVisitor, bit_targets::QasmExpression) = evaluate_bits(v, [bit_targets])


function (v::Qasm2CliffordVisitor)(program_expr::QasmExpression)
    if head(program_expr) == :program
        for expr in program_expr.args
            head(expr) == :end && return
            v(expr)
        end
    elseif head(program_expr) == :version
        v(program_expr.args[1]) == 3.0 || throw(QasmVisitorError("only OpenQASM 3.0 is supported."))
    elseif head(program_expr) == :qubit_declaration
        qubit_name = name(program_expr)
        qubit_size = v(program_expr.args[2])
        qubit_defs(v)[qubit_name] = Qubit(qubit_name, qubit_size)
        qubit_mapping(v)[qubit_name] = collect(qubit_count(v) : qubit_count(v) + qubit_size - 1)
        for qubit_i in 0:qubit_size-1
            qubit_mapping(v)["$qubit_name[$qubit_i]"] = [qubit_count(v) + qubit_i]
        end
        v.qubit_count += qubit_size
    elseif head(program_expr) == :classical_declaration
        init, var_type = declaration_init(v, program_expr)
        var_type isa SizedBitVector || throw(QasmVisitorError("unsupported classical variable type: $var_type."))
        var_name = if head(program_expr.args[2]) == :identifier
            name(program_expr.args[2])
        elseif head(program_expr.args[2]) == :classical_assignment
            name(program_expr.args[2].args[1].args[2])
        else
            throw(QasmVisitorError("invalid classical variable declaration: $program_expr."))
        end
        bit_size = max(v(var_type.size), 1)
        bit_mapping(v)[var_name] = collect(bit_count(v) : bit_count(v) + bit_size - 1)
        for bit_i in 0:bit_size-1
            bit_mapping(v)["$var_name[$bit_i]"] = [bit_count(v) + bit_i]
        end
        v.bit_count += bit_size
        head(program_expr.args[2]) == :classical_assignment && v(program_expr.args[2])
    elseif head(program_expr) == :gate_call
        gate_name = name(program_expr)
        hasgate(v, gate_name) || throw(QasmVisitorError("gate $gate_name not defined!"))
        gate_def = gate_defs(v)[gate_name]
        gate_targets, _, _, _ = process_gate_targets(v, program_expr, gate_def)
        longest, gate_targets = splat_gate_targets(gate_targets)
        for splatted_ix in 1:longest
            args = [gate_targets[i][splatted_ix]+1 for i in eachindex(gate_targets)]
            push!(v, gate_def.body(args...))
        end
    elseif head(program_expr) == :classical_assignment
        left_hand_side  = program_expr.args[1].args[2]
        right_hand_side = program_expr.args[1].args[3]
        head(right_hand_side) == :measure || throw(QasmVisitorError("unsupported assignment operation: $program_expr."))
        bits = evaluate_bits(v, left_hand_side)
        qubits = evaluate_qubits(v, right_hand_side.args[1])
        length(qubits) == length(bits) || throw(QasmVisitorError("measurement source and destination sizes must match."))
        for (q, c) in zip(qubits, bits)
            push!(v, sMZ(q+1, c+1))
        end
    elseif head(program_expr) == :measure
        qubits_to_measure = evaluate_qubits(v, program_expr.args[1])
        for q in qubits_to_measure
            push!(v, sMZ(q+1))
        end
    elseif head(program_expr) == :reset
        qubits_to_reset = evaluate_qubits(v, program_expr.args[1])
        for q in qubits_to_reset
            push!(v, sMRZ(q+1))
        end
    elseif head(program_expr) == :unary_op
        op  = program_expr.args[1]
        arg = v(program_expr.args[2])
        return evaluate_unary_op(op, arg)
    elseif head(program_expr) ∈ (:integer_literal, :float_literal, :string_literal, :complex_literal, :irrational_literal, :boolean_literal, :duration_literal)
        return program_expr.args[1]
    elseif head(program_expr) == :array_literal
        return map(v, program_expr.args)
    elseif head(program_expr) == :range
        return StepRange(v(program_expr.args)...)
    else
        throw(QasmVisitorError("cannot visit expression $program_expr."))
    end
    return v
end