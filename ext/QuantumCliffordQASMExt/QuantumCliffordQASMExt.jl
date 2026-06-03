module QuantumCliffordQASMExt

import QuantumClifford: parse_qasm3, read_qasm3
include("CliffordQasmVisitor.jl")

parse_qasm3(qasm::AbstractString) = Qasm2CliffordVisitor()(parse_qasm(qasm)).instructions
read_qasm3(path::AbstractString) = parse_qasm3(read(path, String))

end