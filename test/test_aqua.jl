using Aqua
using QuantumClifford
using Nemo

Aqua.test_all(QuantumClifford, ambiguities=false)
Aqua.test_ambiguities(QuantumClifford, broken=true)