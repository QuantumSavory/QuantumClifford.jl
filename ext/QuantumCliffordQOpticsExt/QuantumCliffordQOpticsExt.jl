module QuantumCliffordQOpticsExt

using QuantumClifford
import QuantumClifford: mul_left!, mul_right!
using QuantumOpticsBase
using Graphs
using DocStringExtensions
import QuantumOpticsBase: Ket, Operator

const _b2 = SpinBasis(1//2)
const _l0 = spinup(_b2)
const _l1 = spindown(_b2)
const _s₊ = (_l0+_l1)/√2
const _s₋ = (_l0-_l1)/√2
const _i₊ = (_l0+im*_l1)/√2
const _i₋ = (_l0-im*_l1)/√2
const _σ₊ = sigmap(_b2)
const _σ₋ = sigmam(_b2)
const _l00 = projector(_l0)
const _l11 = projector(_l1)
const _id = identityoperator(_b2)
const _z = sigmaz(_b2)
const _x = sigmax(_b2)
const _y = sigmay(_b2)
const _Id = identityoperator(_b2)
const _hadamard = (sigmaz(_b2)+sigmax(_b2))/√2
const _cnot = _l00⊗_Id + _l11⊗_x
const _cphase = _l00⊗_Id + _l11⊗_z
const _phase = _l00 + im*_l11
const _iphase = _l00 - im*_l11

function stab_to_ket(s)
    r,c = size(stabilizerview(s))
    r==c || throw(ArgumentError("""
        There was an attempt to convert a non-square tableaux into a `QuantumOptics.Ket`.
        The Stabilizer tableau has to be square (i.e. pure), otherwise the conversion is impossible.
    """))
    graph, hadamard_idx, iphase_idx, flips_idx = graphstate(s)
    ket = tensor(fill(copy(_s₊),c)...) # TODO fix this is UGLY
    for (;src,dst) in edges(graph)
        apply!(ket, [src,dst], _cphase)
    end
    for i in flips_idx
        apply!(ket, [i], _z)
    end
    for i in iphase_idx
        apply!(ket, [i], _phase)
    end
    for i in hadamard_idx
        apply!(ket, [i], _hadamard)
    end
    ket
end

"""
$TYPEDSIGNATURES

Convert a stabilizer state to a ket representation.

Requires that `QuantumOptics` is loaded.

```jldoctest
julia> using QuantumClifford, QuantumOpticsBase

julia> Ket(S"XX ZZ")
Ket(dim=4)
  basis: [Spin(1/2) ⊗ Spin(1/2)]
 0.7071067811865474 + 0.0im
                0.0 + 0.0im
                0.0 + 0.0im
 0.7071067811865474 + 0.0im
```"""
Ket(s::QuantumClifford.AbstractStabilizer) = stab_to_ket(s)

function stabmix_to_densityop(s::StabMixture)
    ρ₀ = zero(dm(Ket(s.stab)))
    for ((Pₗᵇⁱᵗˢ,Pᵣᵇⁱᵗˢ), χ) in s.destabweights
        ρ̃ = dm(Ket(s.stab))
        destab = destabilizerview(s.stab)
        p = zero(destab[1])
        for (i,b) in enumerate(Pₗᵇⁱᵗˢ)
            b && mul_left!(p, destab, i)
        end
        ρ̃ = Operator(p)*ρ̃
        p = zero(destab[1])
        for (i,b) in enumerate(Pᵣᵇⁱᵗˢ)
            b && mul_left!(p, destab, i)
        end
        ρ̃ = ρ̃*Operator(p)
        ρ₀ += χ*ρ̃
    end
    return ρ₀
end

"""
$TYPEDSIGNATURES

"""
Operator(s::StabMixture) = stabmix_to_densityop(s)


"""
$TYPEDSIGNATURES

Convert a `QuantumClifford.PauliOperator` to `QuantumOptics.Operator`.

```jldoctest
julia> Operator(P"Y") |> dense
Operator(dim=2x2)
  basis: Spin(1/2)
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im

julia> Operator(P"-I") |> dense
Operator(dim=2x2)
  basis: Spin(1/2)
 -1.0+0.0im   0.0+0.0im
  0.0+0.0im  -1.0+0.0im
```
"""
function Operator(p::PauliOperator)
    toop = Dict((false,false)=>_Id,(true,false)=>_x,(false,true)=>_z,(true,true)=>_y)
    return (im)^p.phase[] * ⊗((toop[p[i]] for i in eachindex(p))...)
end

"""
$TYPEDSIGNATURES

Convert a `QuantumClifford.UnitaryPauliChannel` to `QuantumOptics.Operator`.

```jldoctest
julia> pcT
A unitary Pauli channel P = ∑ ϕᵢ Pᵢ with the following branches:
with ϕᵢ | Pᵢ
 0.853553+0.353553im | + _
 0.146447-0.353553im | + Z

julia> Operator(pcT) |> dense
Operator(dim=2x2)
  basis: Spin(1/2)
 1.0+0.0im       0.0+0.0im
 0.0+0.0im  0.707107+0.707107im
```
"""
function Operator(p::UnitaryPauliChannel)
    return sum(w*Operator(p) for (p,w) in zip(p.paulis,p.weights))
end


# TODO you need to decide on big or small endian -- what would make it match QuantumOpticsBase?
function bitstring_to_stabilizer(bitstring::Integer, n::Int)
    s = one(Stabilizer,n)
    for i in 1:n
        if bitstring & (1<<(n-i)) != 0
            s.tab.phases[i] = 0x2
        end
    end
    s
end

# TODO this is not yet functional as it misses the phase
function cliff_to_unitary(cliff)
    basis = [bitstring_to_stabilizer(b,nqubits(cliff)) for b in 0:2^nqubits(cliff)-1]
    sum((projector(stab_to_ket(b),stab_to_ket(cliff*b)') for b in basis))
end

end
