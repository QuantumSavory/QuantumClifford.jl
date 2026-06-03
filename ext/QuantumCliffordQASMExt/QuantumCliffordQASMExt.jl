module QuantumCliffordQASMExt

using Quasar
using QuantumClifford
include("CliffordQasmVisitor.jl")

parse_qasm3(qasm::String) = Qasm2CliffordVisitor()(parse_qasm(qasm)).instructions
read_qasm3(filename::String) = parse_qasm3(read(filename, String))

end