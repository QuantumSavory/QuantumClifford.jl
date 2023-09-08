using LinearAlgebra: Adjoint
using QuantumClifford: Tableau

"""
TableauCUDA is the type of Tableau with its data stored in CUDA arrays
TableauAdj is the type of Tableau with its data stored in an adjoint of a CUDA arrays

rest of the fastrow, fastcolumn functions are implemented in QuantumClifford.jl
"""

# todo when we use to_gpu, to_cpu the efffect of fastrow, fastcolumn disappears
fastrow(t::TableauCUDA) = t
fastrow(t::TableauAdj) = Tableau(t.phases, t.nqubits, CuArray(t.xzs))
fastcolumn(t::TableauCUDA) = Tableau(t.phases, t.nqubits, CuArray(t.xzs')')
fastcolumn(t::TableauAdj) = t