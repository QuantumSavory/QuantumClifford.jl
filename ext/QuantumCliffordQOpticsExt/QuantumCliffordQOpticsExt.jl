module QuantumCliffordQOpticsExt

using QuantumClifford
import QuantumClifford: mul_left!, mul_right!, dot
using QuantumOpticsBase
using Graphs
using DocStringExtensions
import QuantumOpticsBase: Ket, Operator
using LinearAlgebra

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

"""
$TYPEDSIGNATURES

Convert a stabilizer tableau to a density matrix corresponding to the given state.
"""
Operator(s::QuantumClifford.AbstractStabilizer) = dm(Ket(s)) # TODO support mixed stabilizer states

function genstab_to_densityop(s::GeneralizedStabilizer)
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
Operator(s::GeneralizedStabilizer) = genstab_to_densityop(s)


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

function cliff_to_unitary(cliff)
    n = nqubits(cliff)
    b = bell(n, localorder=true)
    apply!(b, cliff, 1:n)
    ψ = Ket(b)
    Operator(SpinBasis(1//2)^n,reshape(ψ.data * sqrt(2)^n, (2^n,2^n)))
end

"""
$TYPEDSIGNATURES

Convert a `QuantumClifford.CliffordOperator` to `QuantumOptics.Operator`.

```jldoctest
julia> Operator(tHadamard)
Operator(dim=2x2)
  basis: Spin(1/2)
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im

julia> Operator(tId1⊗tHadamard)
Operator(dim=4x4)
  basis: [Spin(1/2) ⊗ Spin(1/2)]
 0.707107+0.0im       0.0+0.0im   0.707107+0.0im        0.0+0.0im
      0.0+0.0im  0.707107+0.0im        0.0+0.0im   0.707107+0.0im
 0.707107+0.0im       0.0+0.0im  -0.707107+0.0im        0.0+0.0im
      0.0+0.0im  0.707107+0.0im        0.0+0.0im  -0.707107+0.0im

julia> Operator(tCNOT)
Operator(dim=4x4)
  basis: [Spin(1/2) ⊗ Spin(1/2)]
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
```

This conversion expects a dense tableau of type [`CliffordOperator`](@ref) as input.
If you are working with some of the implicit (a.k.a. small/sparse/symbolic) types of gates
(e.g. [`sCNOT`](@ref)), you need to first convert them to `CliffordOperator`.

```jldoctest
julia> Operator(CliffordOperator(sHadamard))
Operator(dim=2x2)
  basis: Spin(1/2)
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im

julia> Operator(CliffordOperator(sHadamard(1), 1))
Operator(dim=2x2)
  basis: Spin(1/2)
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im

julia> Operator(CliffordOperator(sHadamard(1), 2))
Operator(dim=4x4)
  basis: [Spin(1/2) ⊗ Spin(1/2)]
 0.707107+0.0im   0.707107+0.0im       0.0+0.0im        0.0+0.0im
 0.707107+0.0im  -0.707107+0.0im       0.0+0.0im        0.0+0.0im
      0.0+0.0im        0.0+0.0im  0.707107+0.0im   0.707107+0.0im
      0.0+0.0im        0.0+0.0im  0.707107+0.0im  -0.707107+0.0im

julia> Operator(CliffordOperator(sHadamard(2), 2))
Operator(dim=4x4)
  basis: [Spin(1/2) ⊗ Spin(1/2)]
 0.707107+0.0im       0.0+0.0im   0.707107+0.0im        0.0+0.0im
      0.0+0.0im  0.707107+0.0im        0.0+0.0im   0.707107+0.0im
 0.707107+0.0im       0.0+0.0im  -0.707107+0.0im        0.0+0.0im
      0.0+0.0im  0.707107+0.0im        0.0+0.0im  -0.707107+0.0im

julia> Operator(CliffordOperator(sHadamard(2), 1, compact=true))
Operator(dim=2x2)
  basis: Spin(1/2)
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im
```
"""
function Operator(c::CliffordOperator)
    cliff_to_unitary(c)
end

"""
The inner product of two [`GeneralizedStabilizer`](@ref) states, `sm₁` and `sm₂`.

```jldoctest
julia> using QuantumOpticsBase; using LinearAlgebra; # hide

julia> sm = GeneralizedStabilizer(S"X")
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
+ X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 1.0+0.0im | + _ | + _

julia> apply!(sm, pcT)
A mixture ∑ ϕᵢⱼ Pᵢ ρ Pⱼ† where ρ is
𝒟ℯ𝓈𝓉𝒶𝒷
+ Z
𝒮𝓉𝒶𝒷
+ X
with ϕᵢⱼ | Pᵢ | Pⱼ:
 0.0+0.353553im | + _ | + Z
 0.0-0.353553im | + Z | + _
 0.853553+0.0im | + _ | + _
 0.146447+0.0im | + Z | + Z

julia> dot(sm, sm)
0.9999999999999994
```
"""
function LinearAlgebra.dot(sm₁::GeneralizedStabilizer, sm₂::GeneralizedStabilizer)
    return real(tr(Operator(sm₁)' * Operator(sm₂)))
end

end

##
