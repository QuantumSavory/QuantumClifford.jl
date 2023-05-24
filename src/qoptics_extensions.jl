##
# Extensions
##

"""Convert a stabilizer state to a ket representation.

Requires that `QuantumOptics` is loaded.

```jldoctest
julia> using QuantumClifford, QuantumOpticsBase

julia> stab_to_ket(S"XX ZZ")
Ket(dim=4)
  basis: [Spin(1/2) ⊗ Spin(1/2)]
 0.7071067811865474 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865474 + 0.0im
```"""
function stab_to_ket end

function cliff_to_unitary end
