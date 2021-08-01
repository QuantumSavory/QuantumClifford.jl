using BenchmarkTools
using QuantumClifford

const SUITE = BenchmarkGroup()

SUITE["pauli"] = BenchmarkGroup(["pauli"])
SUITE["pauli"]["mul"] = BenchmarkGroup(["multiplication"])
SUITE["pauli"]["mul"]["100"] = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=random_pauli(100);b=random_pauli(100))