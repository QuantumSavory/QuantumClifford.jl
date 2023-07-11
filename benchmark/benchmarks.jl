using BenchmarkTools
using Pkg
using StableRNGs
using QuantumClifford
using Nemo

const SUITE = BenchmarkGroup()

rng = StableRNG(42)

M = Pkg.Operations.Context().env.manifest
V = M[findfirst(v -> v.name == "QuantumClifford", M)].version

# Due to breaking changes in v0.5.0
const tableauCNOT = C"XX
                      IX
                      ZI
                      ZZ"

const s500 = random_stabilizer(rng,500)
const md500 = MixedDestabilizer(random_destabilizer(rng,500),250)
const p500 = random_pauli(rng,500)
const c500 = random_clifford(rng,500)
const cnots250 = tensor_pow(tableauCNOT,250)
const pa1 = random_pauli(rng,100)
const pb1 = random_pauli(rng,100)
const pa2 = random_pauli(rng,1000)
const pb2 = random_pauli(rng,1000)
const pa3 = random_pauli(rng,100_000)
const pb3 = random_pauli(rng,100_000)
const pa4 = random_pauli(rng,20_000_000)
const pb4 = random_pauli(rng,20_000_000)

SUITE["pauli"] = BenchmarkGroup(["pauli"])
# Multiplication of multiqubit Pauli operators of various length.
SUITE["pauli"]["mul"] = BenchmarkGroup(["multiplication"])
SUITE["pauli"]["mul"]["100"]    = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=copy(pa1);b=pb1)
SUITE["pauli"]["mul"]["1000"]   = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=copy(pa2);b=pb2)
SUITE["pauli"]["mul"]["100000"] = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=copy(pa3);b=pb3)
SUITE["pauli"]["mul"]["20000000"] = @benchmarkable QuantumClifford.mul_left!(a,b) setup=(a=copy(pa4);b=pb4)

SUITE["stabilizer"] = BenchmarkGroup(["stabilizer"])
# Canonicalization of random or diagonal tableaux. All three types of canonicalization are benchmarked.
SUITE["stabilizer"]["canon"] = BenchmarkGroup(["canonicalization"])
SUITE["stabilizer"]["canon"]["cano500"] = @benchmarkable canonicalize!(s) setup=(s=copy(s500)) evals=1
SUITE["stabilizer"]["canon"]["gott500"] = @benchmarkable canonicalize_gott!(s) setup=(s=copy(s500)) evals=1
SUITE["stabilizer"]["canon"]["rref500"] = @benchmarkable canonicalize_rref!(s) setup=(s=copy(s500)) evals=1
SUITE["stabilizer"]["canon"]["md_cano500"] = @benchmarkable canonicalize!(md) setup=(md=copy(md500)) evals=1
SUITE["stabilizer"]["canon"]["md_rref500"] = @benchmarkable canonicalize_rref!(md) setup=(md=copy(md500)) evals=1
SUITE["stabilizer"]["canon"]["diag_cano500"] = @benchmarkable canonicalize!(s) setup=(s=one(Stabilizer,500))  evals=1
SUITE["stabilizer"]["canon"]["diag_gott500"] = @benchmarkable canonicalize_gott!(s) setup=(s=one(Stabilizer,500)) evals=1
SUITE["stabilizer"]["canon"]["diag_rref500"] = @benchmarkable canonicalize_rref!(s) setup=(s=one(Stabilizer,500)) evals=1
# Tensor products of tableaux.
SUITE["stabilizer"]["tensor"] = BenchmarkGroup(["tensor product"])
SUITE["stabilizer"]["tensor"]["pow5_20"] = @benchmarkable tensor_pow(s,20) setup=(s=random_stabilizer(5))
SUITE["stabilizer"]["tensor"]["diag_pow5_20"] = @benchmarkable tensor_pow(s,20) setup=(s=one(Stabilizer,500))
# Projection operations.
SUITE["stabilizer"]["project"] = BenchmarkGroup(["project", "measure"])
SUITE["stabilizer"]["project"]["stabilizer"]   = @benchmarkable project!(s,p) setup=(s=copy(s500);p=copy(p500)) evals=1
SUITE["stabilizer"]["project"]["destabilizer"] = @benchmarkable project!(s,p) setup=(s=copy(md500);p=copy(p500)) evals=1
# Partial traces.
SUITE["stabilizer"]["trace"] = BenchmarkGroup(["partial trace", "erase"])
SUITE["stabilizer"]["trace"]["stabilizer"]   = @benchmarkable traceout!(s,[2,4,70]) setup=(s=copy(s500);p=copy(p500)) evals=1
SUITE["stabilizer"]["trace"]["destabilizer"] = @benchmarkable traceout!(s,[2,4,70]) setup=(s=copy(md500);p=copy(p500)) evals=1

SUITE["clifford"] = BenchmarkGroup(["clifford"])
# Tensor products of Clifford operators # TODO finish functionality and add benchmarks
#SUITE["clifford"]["tensor"] = BenchmarkGroup(["tensor product"])
# Applications of various Clifford operators to tableaux.
SUITE["clifford"]["dense"] = BenchmarkGroup(["dense"])
# Both the operator and the state are random dense matrices.
SUITE["clifford"]["dense"]["dense500_on_dense500_stab"]   = @benchmarkable apply!(s,c) setup=(s=copy(s500);c=c500) evals=1
SUITE["clifford"]["dense"]["dense500_on_dense500_destab"] = @benchmarkable apply!(s,c) setup=(s=copy(md500);c=c500) evals=1
# The operator is a small 2-qubit CNOT operator.
SUITE["clifford"]["dense"]["cnot_on_dense500_stab"]   = @benchmarkable apply!(s,c,[30,200]) setup=(s=copy(s500);c=tableauCNOT) evals=1
SUITE["clifford"]["dense"]["cnot_on_dense500_destab"] = @benchmarkable apply!(s,c,[30,200]) setup=(s=copy(md500);c=tableauCNOT) evals=1
# The state is a diagonal matrix (and the operator is dense large matrix or just a dense 2-qubit CNOT matrix).
SUITE["clifford"]["dense"]["dense500_on_diag500_stab"]   = @benchmarkable apply!(s,c) setup=(s=one(Stabilizer,500);c=c500) evals=1
SUITE["clifford"]["dense"]["dense500_on_diag500_destab"] = @benchmarkable apply!(s,c) setup=(s=one(Destabilizer,500);c=c500) evals=1
SUITE["clifford"]["dense"]["cnot_on_diag500_stab"]   = @benchmarkable apply!(s,c,[30,200]) setup=(s=one(Stabilizer,500);c=tableauCNOT) evals=1
SUITE["clifford"]["dense"]["cnot_on_diag500_destab"] = @benchmarkable apply!(s,c,[30,200]) setup=(s=MixedDestabilizer(one(Destabilizer,500),250);c=tableauCNOT) evals=1
# The operator is a mostly diagonal matrix (250 CNOTs), while the state is a random dense matrix.
SUITE["clifford"]["dense"]["cnot250_on_dense500_stab"]   = @benchmarkable apply!(s,c) setup=(s=copy(s500);c=cnots250) evals=1
SUITE["clifford"]["dense"]["cnot250_on_dense500_destab"] = @benchmarkable apply!(s,c) setup=(s=copy(md500);c=cnots250) evals=1
# The operator is mostly diagonal and the state is diagonal.
SUITE["clifford"]["dense"]["cnot250_on_diag500_stab"]   = @benchmarkable apply!(s,c) setup=(s=one(Stabilizer,500);c=cnots250) evals=1
SUITE["clifford"]["dense"]["cnot250_on_diag500_destab"] = @benchmarkable apply!(s,c) setup=(s=MixedDestabilizer(one(Destabilizer,500),250);c=cnots250) evals=1
# Tensor products of symbolic Clifford operators # TODO finish functionality and add benchmarks
SUITE["clifford"]["symbolic"] = BenchmarkGroup(["symbolic"])
SUITE["clifford"]["symbolic"]["cnot_on_dense500_stab"]   = @benchmarkable apply!(s,c) setup=(s=copy(s500);c=sCNOT(30,200)) evals=1
SUITE["clifford"]["symbolic"]["cnot_on_dense500_destab"] = @benchmarkable apply!(s,c) setup=(s=copy(md500);c=sCNOT(30,200)) evals=1
SUITE["clifford"]["symbolic"]["cnot_on_diag500_stab"]   = @benchmarkable apply!(s,c) setup=(s=one(Stabilizer,500);c=sCNOT(30,200)) evals=1
SUITE["clifford"]["symbolic"]["cnot_on_diag500_destab"] = @benchmarkable apply!(s,c) setup=(s=MixedDestabilizer(one(Destabilizer,500),250);c=sCNOT(30,200)) evals=1
SUITE["clifford"]["symbolic"]["cnot250_on_dense500_stab"]   = @benchmarkable for i in 1:250 apply!(s,sCNOT(2i-1,2i)) end setup=(s=copy(s500)) evals=1
SUITE["clifford"]["symbolic"]["cnot250_on_dense500_destab"] = @benchmarkable for i in 1:250 apply!(s,sCNOT(2i-1,2i)) end setup=(s=copy(md500)) evals=1
SUITE["clifford"]["symbolic"]["cnot250_on_diag500_stab"]   = @benchmarkable for i in 1:250 apply!(s,sCNOT(2i-1,2i)) end setup=(s=one(Stabilizer,500)) evals=1
SUITE["clifford"]["symbolic"]["cnot250_on_diag500_destab"] = @benchmarkable for i in 1:250 apply!(s,sCNOT(2i-1,2i)) end setup=(s=MixedDestabilizer(one(Destabilizer,500),250)) evals=1


if V > v"0.7.0"

SUITE["circuitsim"] = BenchmarkGroup(["circuitsim"])
function x_diag_circuit(csize)
    circuit = []
    for i in 1:csize
        push!(circuit, sHadamard(i))
        push!(circuit, sCNOT(i, csize+1))
    end
    push!(circuit, sMZ(csize+1,1))
    return circuit
end
SUITE["circuitsim"]["pftrajectories"] = BenchmarkGroup(["pftrajectories"])
SUITE["circuitsim"]["pftrajectories"]["q101_r1"]      = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(   1,  101, 1); circuit=x_diag_circuit( 100)) evals=1
SUITE["circuitsim"]["pftrajectories"]["q1001_r1"]     = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(   1, 1001, 1); circuit=x_diag_circuit(1000)) evals=1
SUITE["circuitsim"]["pftrajectories"]["q1001_r100"]   = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame( 100, 1001, 1); circuit=x_diag_circuit(1000)) evals=1
SUITE["circuitsim"]["pftrajectories"]["q1001_r10000"] = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(1000, 1001, 1); circuit=x_diag_circuit(1000)) evals=1
SUITE["circuitsim"]["mctrajectories"] = BenchmarkGroup(["mctrajectories"])
SUITE["circuitsim"]["mctrajectories"]["q101_r1"]    = @benchmarkable mctrajectory!(state, circuit) setup=(state=Register(one(Stabilizer,  101), [false]); circuit=x_diag_circuit( 100)) evals=1
SUITE["circuitsim"]["mctrajectories"]["q1001_r1"]   = @benchmarkable mctrajectory!(state, circuit) setup=(state=Register(one(Stabilizer, 1001), [false]); circuit=x_diag_circuit(1000)) evals=1

SUITE["circuitsim"]["pftrajectories_union"] = BenchmarkGroup(["pftrajectories_union"])
veccirc = Vector{Union{sHadamard,sCNOT,sMZ}}
SUITE["circuitsim"]["pftrajectories_union"]["q101_r1"]      = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(   1,  101, 1); circuit=veccirc(x_diag_circuit( 100))) evals=1
SUITE["circuitsim"]["pftrajectories_union"]["q1001_r1"]     = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(   1, 1001, 1); circuit=veccirc(x_diag_circuit(1000))) evals=1
SUITE["circuitsim"]["pftrajectories_union"]["q1001_r100"]   = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame( 100, 1001, 1); circuit=veccirc(x_diag_circuit(1000))) evals=1
SUITE["circuitsim"]["pftrajectories_union"]["q1001_r10000"] = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(1000, 1001, 1); circuit=veccirc(x_diag_circuit(1000))) evals=1
SUITE["circuitsim"]["mctrajectories_union"] = BenchmarkGroup(["mctrajectories_union"])
SUITE["circuitsim"]["mctrajectories_union"]["q101_r1"]    = @benchmarkable mctrajectory!(state, circuit) setup=(state=Register(one(Stabilizer,  101), [false]); circuit=veccirc(x_diag_circuit( 100))) evals=1
SUITE["circuitsim"]["mctrajectories_union"]["q1001_r1"]   = @benchmarkable mctrajectory!(state, circuit) setup=(state=Register(one(Stabilizer, 1001), [false]); circuit=veccirc(x_diag_circuit(1000))) evals=1

SUITE["circuitsim"]["pftrajectories_sumtype"] = BenchmarkGroup(["pftrajectories_sumtype"])
SUITE["circuitsim"]["pftrajectories_sumtype"]["q101_r1"]      = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(   1,  101, 1); circuit=compactify_circuit(x_diag_circuit( 100))) evals=1
SUITE["circuitsim"]["pftrajectories_sumtype"]["q1001_r1"]     = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(   1, 1001, 1); circuit=compactify_circuit(x_diag_circuit(1000))) evals=1
SUITE["circuitsim"]["pftrajectories_sumtype"]["q1001_r100"]   = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame( 100, 1001, 1); circuit=compactify_circuit(x_diag_circuit(1000))) evals=1
SUITE["circuitsim"]["pftrajectories_sumtype"]["q1001_r10000"] = @benchmarkable pftrajectories(state,circuit) setup=(state=PauliFrame(1000, 1001, 1); circuit=compactify_circuit(x_diag_circuit(1000))) evals=1
SUITE["circuitsim"]["pftrajectories_sumtype"]["q1001_r10000_fastrow"] = @benchmarkable pftrajectories(state,circuit) setup=(state=fastrow(PauliFrame(1000, 1001, 1)); circuit=compactify_circuit(x_diag_circuit(1000))) evals=1
SUITE["circuitsim"]["mctrajectories_sumtype"] = BenchmarkGroup(["mctrajectories_sumtype"])
SUITE["circuitsim"]["mctrajectories_sumtype"]["q101_r1"]    = @benchmarkable mctrajectory!(state, circuit) setup=(state=Register(one(Stabilizer,  101), [false]); circuit=compactify_circuit(x_diag_circuit( 100))) evals=1
SUITE["circuitsim"]["mctrajectories_sumtype"]["q1001_r1"]   = @benchmarkable mctrajectory!(state, circuit) setup=(state=Register(one(Stabilizer, 1001), [false]); circuit=compactify_circuit(x_diag_circuit(1000))) evals=1

end
