# OpenQASM 3 Import

QuantumClifford can import circuits described in
[OpenQASM 3](https://openqasm.com/versions/3.0/language/) via the optional
[Quasar.jl](https://github.com/kshyatt-aws/Quasar.jl) extension.

## Quick start

First install Quasar.jl:

```julia-repl
pkg> add Quasar
```

Then load both packages — loading `Quasar` automatically activates the extension:

```julia
using QuantumClifford
using Quasar                     # activates QuantumCliffordQuasarExt

circuit = read_qasm3("bell.qasm")
frames  = pftrajectories(circuit; trajectories=10_000)
```

You can also parse directly from a string:

```julia
src = """
OPENQASM 3.0;

qubit[2] q;
bit[2] c;

h  q[0];
cx q[0], q[1];
c[0] = measure q[0];
c[1] = measure q[1];
"""

circuit = parse_qasm3(src)
# => [sHadamard(1), sCNOT(1, 2), sMZ(1, 1), sMZ(2, 2)]
```

## Index convention

OpenQASM 3 uses **0-based** qubit and bit indices; QuantumClifford uses
**1-based** indices.  The conversion is performed automatically:

| OpenQASM 3       | QuantumClifford |
|------------------|-----------------|
| `q[0]`           | qubit `1`       |
| `q[1]`           | qubit `2`       |
| `c[0]`           | classical bit `1` |

## Supported operations

### Single-qubit Clifford gates

| OpenQASM 3 name | QuantumClifford operation |
|-----------------|--------------------------|
| `id` / `i`      | `sId1(q)`                |
| `x`             | `sX(q)`                  |
| `y`             | `sY(q)`                  |
| `z`             | `sZ(q)`                  |
| `h`             | `sHadamard(q)`           |
| `s`             | `sPhase(q)`              |
| `sdg`           | `sInvPhase(q)`           |

### Two-qubit Clifford gates

| OpenQASM 3 name   | QuantumClifford operation |
|-------------------|--------------------------|
| `cx` / `cnot`     | `sCNOT(q1, q2)`          |
| `cz`              | `sCPHASE(q1, q2)`        |
| `swap`            | `sSWAP(q1, q2)`          |

### Classical operations

| OpenQASM 3 statement      | QuantumClifford operation |
|---------------------------|--------------------------|
| `c[i] = measure q[j]`    | `sMZ(q, c)` — measure qubit `q` in the Z basis and store result in classical bit `c` |
| `reset q[i]`              | `sMRZ(q, 0)` — measure qubit `q` and reset to \|0⟩ |

### Register declarations

Both `qubit[n] name` and `bit[n] name` declarations are supported.
Multiple registers of each kind are allowed; the qubit (and bit) indices are
allocated sequentially in declaration order:

```qasm
OPENQASM 3.0;
qubit[2] a;    // a[0] → qubit 1, a[1] → qubit 2
qubit[2] b;    // b[0] → qubit 3, b[1] → qubit 4
bit[2]   ca;   // ca[0] → bit 1,  ca[1] → bit 2
bit[2]   cb;   // cb[0] → bit 3,  cb[1] → bit 4
```

## Unsupported operations

Any gate or construct that is not in the table above will raise an
`ArgumentError` that includes the gate name and, when available, its source
location.  This includes:

- Non-Clifford gates (`t`, `tdg`, `rx`, `ry`, `rz`, `u`, `u1`, `u2`, `u3`, …)
- Parameterized gate calls with non-zero angle parameters
- Custom gate definitions
- Control flow (`if`, `for`, `while`)
- Classical arithmetic and feed-forward
- Timing, calibration, and pulse-level constructs

Example error message:

```
ArgumentError: Unsupported single-qubit gate 't' at line 5, column 1.
Supported gates: ["h", "id", "s", "sdg", "x", "y", "z"].
```

## API reference

```@docs
parse_qasm3
read_qasm3
```
