using BenchmarkTools
using QuantumClifford
using Nemo

const SUITE = BenchmarkGroup()

SUITE["pauli"] = BenchmarkGroup(["pauli"])
# Multiplication of multiqubit Pauli operators of various length.
SUITE["pauli"]["mul"] = BenchmarkGroup(["multiplication"])
SUITE["pauli"]["mul"]["100"]    = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=random_pauli(100);b=random_pauli(100))
SUITE["pauli"]["mul"]["1000"]   = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=random_pauli(1_000);b=random_pauli(1_000))
SUITE["pauli"]["mul"]["100000"] = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=random_pauli(100_000);b=random_pauli(100_000))

SUITE["stabilizer"] = BenchmarkGroup(["stabilizer"])
# Canonicalization of random or diagonal tableaux. All three types of canonicalization are benchmarked.
SUITE["stabilizer"]["canon"] = BenchmarkGroup(["canonicalization"])
SUITE["stabilizer"]["canon"]["cano500"] = @benchmarkable canonicalize!(s) setup=(s=random_stabilizer(500)) evals=1
SUITE["stabilizer"]["canon"]["gott500"] = @benchmarkable canonicalize_gott!(s) setup=(s=random_stabilizer(500)) evals=1
SUITE["stabilizer"]["canon"]["rref500"] = @benchmarkable canonicalize_rref!(s) setup=(s=random_stabilizer(500)) evals=1
SUITE["stabilizer"]["canon"]["diag_cano500"] = @benchmarkable canonicalize!(s) setup=(s=one(Stabilizer,500))  evals=1
SUITE["stabilizer"]["canon"]["diag_gott500"] = @benchmarkable canonicalize_gott!(s) setup=(s=one(Stabilizer,500)) evals=1
SUITE["stabilizer"]["canon"]["diag_rref500"] = @benchmarkable canonicalize_rref!(s) setup=(s=one(Stabilizer,500)) evals=1
# Tensor products of tableaux. # TODO finish functionality and add benchmarks
#SUITE["stabilizer"]["tensor"] = BenchmarkGroup(["tensor product"])
#SUITE["stabilizer"]["tensor"]["pow5_20"] = @benchmarkable tensor_pow(s,20) setup=(s=random_stabilizer(5))
#SUITE["stabilizer"]["tensor"]["diag_pow5_20"] = @benchmarkable tensor_pow(s,20) setup=(s=one(Stabilizer,500))
# Projection operations.
SUITE["stabilizer"]["project"] = BenchmarkGroup(["project", "measure"])
SUITE["stabilizer"]["project"]["stabilizer"]   = @benchmarkable project!(s,p) setup=(s=random_stabilizer(500);p=random_pauli(500)) evals=1
SUITE["stabilizer"]["project"]["destabilizer"] = @benchmarkable project!(s,p) setup=(s=MixedDestabilizer(random_destabilizer(500).tab,500);p=random_pauli(500)) evals=1 # TODO make a random_mixeddestabilizer and a MixedDestabilizer(s::Destabilizer)
# Partial traces.
SUITE["stabilizer"]["trace"] = BenchmarkGroup(["partial trace", "erase"])
SUITE["stabilizer"]["trace"]["stabilizer"]   = @benchmarkable traceout!(s,[2,4,70]) setup=(s=random_stabilizer(500);p=random_pauli(500)) evals=1
SUITE["stabilizer"]["trace"]["destabilizer"] = @benchmarkable traceout!(s,[2,4,70]) setup=(s=MixedDestabilizer(random_destabilizer(500).tab,500);p=random_pauli(500)) evals=1

SUITE["clifford"] = BenchmarkGroup(["clifford"])
# Tensor products of Clifford operators # TODO finish functionality and add benchmarks
#SUITE["clifford"]["tensor"] = BenchmarkGroup(["tensor product"])
# Applications of various Clifford operators to tableaux.
SUITE["clifford"]["dense"] = BenchmarkGroup(["dense"])
# Both the operator and the state are random dense matrices.
SUITE["clifford"]["dense"]["dense500_on_dense500_stab"]   = @benchmarkable apply!(s,c) setup=(s=random_stabilizer(500);c=random_clifford(500))
SUITE["clifford"]["dense"]["dense500_on_dense500_destab"] = @benchmarkable apply!(s,c) setup=(s=MixedDestabilizer(random_destabilizer(500).tab,500);c=random_clifford(500))
# The operator is a small 2-qubit random operator.
SUITE["clifford"]["dense"]["dense2_on_dense500_stab"]   = @benchmarkable apply!(s,c,[30,200]) setup=(s=random_stabilizer(500);c=(CNOT))
SUITE["clifford"]["dense"]["dense2_on_dense500_destab"] = @benchmarkable apply!(s,c,[30,200]) setup=(s=MixedDestabilizer(random_destabilizer(500).tab,500);c=CliffordOperator(CNOT))
# The state is a diagonal matrix (and the operator is dense large matrix or just a dense 2-qubit matrix).
SUITE["clifford"]["dense"]["dense500_on_diag500_stab"]   = @benchmarkable apply!(s,c) setup=(s=one(Stabilizer,500);c=random_clifford(500))
SUITE["clifford"]["dense"]["dense500_on_diag500_destab"] = @benchmarkable apply!(s,c) setup=(s=one(Destabilizer,500);c=random_clifford(500))
SUITE["clifford"]["dense"]["dense2_on_diag500_stab"]   = @benchmarkable apply!(s,c,[30,200]) setup=(s=one(Stabilizer,500);c=(CNOT))
SUITE["clifford"]["dense"]["dense2_on_diag500_destab"] = @benchmarkable apply!(s,c,[30,200]) setup=(s=one(Destabilizer,500);c=CliffordOperator(CNOT))
# The operator is a mostly diagonal matrix (250 CNOTs), while the state is a random dense matrix.
SUITE["clifford"]["dense"]["cnot250_on_dense500_stab"]   = @benchmarkable apply!(s,c) setup=(s=random_stabilizer(500);c=tensor_pow((CNOT),250))
SUITE["clifford"]["dense"]["cnot250_on_dense500_destab"] = @benchmarkable apply!(s,c) setup=(s=MixedDestabilizer(random_destabilizer(500).tab,500);c=tensor_pow(CliffordOperator(CNOT),250))
# The operator is mostly diagonal and the state is diagonal.
SUITE["clifford"]["dense"]["cnot250_on_diag500_stab"]   = @benchmarkable apply!(s,c) setup=(s=one(Stabilizer,500);c=tensor_pow((CNOT),250))
SUITE["clifford"]["dense"]["cnot250_on_diag500_destab"] = @benchmarkable apply!(s,c) setup=(s=one(Destabilizer,500);c=tensor_pow(CliffordOperator(CNOT),250))
# Tensor products of symbolic Clifford operators # TODO finish functionality and add benchmarks
#SUITE["clifford"]["symbolic"] = BenchmarkGroup(["symbolic"])