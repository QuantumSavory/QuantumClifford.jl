"""$TYPEDSIGNATURES
Parse an OpenQASM 3 program and convert it to a `Vector` of QuantumClifford 
symbolic operations.

Currently supported OpenQASM features include:
- `qubit` and `bit` register declarations
- Single-qubit Clifford gates: `id`, `x`, `y`, `z`, `h`, `s`, `sdg`
- Two-qubit Clifford gates: `cx`, `cz`, `swap`
- Computational-basis measurement into classical bits
- `reset` operation
Unsupported OpenQASM features will throw a `QasmVisitorError` when encountered.

# Example
```jldoctest
julia> qasm = \"\"\"
    OPENQASM 3.0;

    qubit[2] q;
    bit[2] c;

    h q[0];
    cx q[0], q[1];
    c[0] = measure q[0];
    c[1] = measure q[1];
\"\"\";

julia> parse_qasm3(qasm)
4-element Vector{QuantumClifford.AbstractOperation}:
 sHadamard(1)
 sCNOT(1,2)
 sMZ(1, 1)
 sMZ(2, 2)
```

# Supported operation mappings
| OpenQASM 3 | QuantumClifford |
|------------|-----------------|
| `id`       | `sId` |
| `x`        | `sX` |
| `y`        | `sY` |
| `z`        | `sZ` |
| `h`        | `sHadamard` |
| `s`        | `sPhase` |
| `sdg`      | `sInvPhase` |
| `cx`       | `sCNOT` |
| `cz`       | `sCZ` |
| `swap`     | `sSWAP` |
| `measure`  | `sMZ` |
| `reset`    | `sMRZ` |
"""
function parse_qasm3 end

"""$TYPEDSIGNATURES
Read an OpenQASM 3 program from `path` and convert it to a `Vector` of 
QuantumClifford symbolic operations.
"""
function read_qasm3 end