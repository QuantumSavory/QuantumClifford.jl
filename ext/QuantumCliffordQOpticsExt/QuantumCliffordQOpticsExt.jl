module QuantumCliffordQOpticsExt

using QuantumClifford
using QuantumOpticsBase
using Graphs

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

function QuantumClifford.stab_to_ket(s::Stabilizer)
    r,c = size(s)
    @assert r==c "The Stabilizer tableau has to be square"
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
function QuantumClifford.cliff_to_unitary(cliff)
    basis = [bitstring_to_stabilizer(b,nqubits(cliff)) for b in 0:2^nqubits(cliff)-1]
    sum((projector(stab_to_ket(b),stab_to_ket(cliff*b)') for b in basis))
end

# TODO this should be upstreamed to QuantumOpticsBase
using QuantumInterface: embed, basis, dm, AbstractSuperOperator
import QuantumInterface: apply!
using QuantumOpticsBase: Ket, Operator

function apply!(state::Ket, indices, operation::Operator)
    op = basis(state)==basis(operation) ? operation : embed(basis(state), indices, operation)
    state.data = (op*state).data
    state
end

function apply!(state::Operator, indices, operation::Operator)
    op = basis(state)==basis(operation) ? operation : embed(basis(state), indices, operation)
    state.data = (op*state*op').data
    state
end

function apply!(state::Ket, indices, operation::T) where {T<:AbstractSuperOperator}
    apply!(dm(state), indices, operation)
end

function apply!(state::Operator, indices, operation::T) where {T<:AbstractSuperOperator}
    op = basis(state)==basis(operation) ? operation : embed(basis(state), indices, operation)
    state.data = (op*state).data
    state
end

end
