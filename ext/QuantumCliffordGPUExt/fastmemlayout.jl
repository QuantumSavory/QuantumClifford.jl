using LinearAlgebra: Adjoint
using QuantumClifford: Tableau

# when we use to_gpu, to_cpu the efffect of fastrow, fastcolumn disappears
fastrow(t::TableauCUDA) = t
fastrow(t::TableauAdj) = Tableau(t.phases, t.nqubits, CuArray(t.xzs))
fastcolumn(t::TableauCUDA) = Tableau(t.phases, t.nqubits, CuArray(t.xzs')')
fastcolumn(t::TableauAdj) = t