using WeakDepHelpers: @declare_method_is_in_extension

const parse_qasm3_docstring = """
    parse_qasm3(qasm::AbstractString)

Parse an OpenQASM 3 program into a `Vector` of QuantumClifford symbolic circuit
operations. This method is provided by the Quasar.jl extension, so load it with
`import Quasar` before calling `parse_qasm3`.

The first implementation supports the Clifford subset:

| OpenQASM 3 operation | QuantumClifford operation |
|:---------------------|:--------------------------|
| `id`                 | [`sId1`](@ref)            |
| `x`                  | [`sX`](@ref)              |
| `y`                  | [`sY`](@ref)              |
| `z`                  | [`sZ`](@ref)              |
| `h`                  | [`sHadamard`](@ref)       |
| `s`                  | [`sPhase`](@ref)          |
| `sdg`                | [`sInvPhase`](@ref)       |
| `cx`                 | [`sCNOT`](@ref)           |
| `cz`                 | [`sCPHASE`](@ref)         |
| `swap`               | [`sSWAP`](@ref)           |
| `measure`            | [`sMZ`](@ref)             |
| `reset`              | [`sMRZ`](@ref)            |

OpenQASM uses zero-based qubit and bit indices; the returned QuantumClifford
operations use one-based indices.

# Example

```jldoctest
julia> import Quasar

julia> circuit = parse_qasm3(\"\"\"
       OPENQASM 3.0;
       qubit[2] q;
       bit[2] c;
       h q[0];
       cx q[0], q[1];
       c[0] = measure q[0];
       c[1] = measure q[1];
       \"\"\")
4-element Vector{QuantumClifford.AbstractOperation}:
 sHadamard(1)
 sCNOT(1,2)
 sMZ(1, 1)
 sMZ(2, 2)
```

See also: [`read_qasm3`](@ref), [`pftrajectories`](@ref).
"""

const read_qasm3_docstring = """
    read_qasm3(path::AbstractString)

Read an OpenQASM 3 program from `path` and convert it to a `Vector` of
QuantumClifford symbolic circuit operations.

This method is provided by the Quasar.jl extension, so load it with
`import Quasar` before calling `read_qasm3`.

See also: [`parse_qasm3`](@ref).
"""

@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS parse_qasm3 (:Quasar,) parse_qasm3_docstring
@declare_method_is_in_extension QuantumClifford.WEAKDEP_METHOD_ERROR_HINTS read_qasm3 (:Quasar,) read_qasm3_docstring
