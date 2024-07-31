const all_single_qubit_patterns = (
    (true, false, false, true), # X, Z ↦ X, Z
    (false, true, true, true),  # X, Z ↦ Z, Y
    (true, true, true, false),  # X, Z ↦ Y, X
    (false, true, true, false), # X, Z ↦ Z, X - Hadamard
    (true, false, true, true),  # X, Z ↦ X, Y
    (true, true, false, true)   # X, Z ↦ Y, Z - Phase
)

"""Generate a symbolic single-qubit gate given its index. Optionally, set non-trivial phases.

```jldoctest
julia> enumerate_single_qubit_gates(6)
sPhase on qubit 1
X₁ ⟼ + Y
Z₁ ⟼ + Z

julia> enumerate_single_qubit_gates(6, qubit=2, phases=(true, true))
SingleQubitOperator on qubit 2
X₁ ⟼ - Y
Z₁ ⟼ - Z
```

See also: [`enumerate_cliffords`](@ref)."""
function enumerate_single_qubit_gates(index; qubit=1, phases::Tuple{Bool,Bool}=(false,false))
    @assert index<=6 "Only 6 single-qubit gates exit, up to the choice of phases"
    if phases==(false,false)
        if index==4
            return sHadamard(qubit)
        elseif index==6
            return sPhase(qubit)
        else
            return SingleQubitOperator(qubit, all_single_qubit_patterns[index]..., false, false)
        end
    else
        if index==1
            if     (phases[1], phases[2]) == (false, false)
                return sId1(qubit)
            elseif (phases[1], phases[2]) == (false,  true)
                return sX(qubit)
            elseif (phases[1], phases[2]) == (true,  false)
                return sZ(qubit)
            else
                return sY(qubit)
            end
        else
            return SingleQubitOperator(qubit, all_single_qubit_patterns[index]..., phases...)
        end
    end
end

"""The size of the Clifford group `𝒞` over a given number of qubits, possibly modulo the phases.

For n qubits, not accounting for phases is `2ᵏΠⱼ₌₁ⁿ(4ʲ-1)` where `k = n²`. There are `4ⁿ` different phase configurations.

```jldoctest
julia> clifford_cardinality(7)
457620995529680351512370381586432000
```

Calculate the size of the Symplectic group `𝐒p(2n) ≡ 𝒞ₙ/𝒫ₙ`, where `𝒫ₙ` is the Pauli group over `n` qubits, by setting `phases = false`.

```jldoctest
julia> clifford_cardinality(7, phases=false)
27930968965434591767112450048000
```

See also: [`enumerate_cliffords`](@ref).
"""
function clifford_cardinality(n::Int; phases=true)
    base = BigInt(2)^(n^2) * prod(1:n) do j; BigInt(4)^j-1 end
    if phases
        base*BigInt(2)^(2n)
    else
        base
    end
end

@inline function _findanticommGS(basis, start, n, padded_n, δn)
    for i in δn+start+1:padded_n
        comm(basis, δn+start, i)==0x1 && return i
    end
    for i in padded_n+δn+start:2padded_n+1
        comm(basis, δn+start, i)==0x1 && return i
    end
    # the following happens only if the input is P"X___..."
    rowswap!(basis, δn+start, 2padded_n+1; phases=Val(false))
    _findanticommGS(basis, start, n, padded_n, δn)
end

@inline function _eliminateGS(basis, start, n, padded_n, δn)
    for i in δn+start+1:padded_n
        x = comm(basis, δn+start, i)==0x1
        z = comm(basis, padded_n+δn+start, i)==0x1
        z && mul_left!(basis, i, δn+start; phases=Val(false))
        x && mul_left!(basis, i, padded_n+δn+start; phases=Val(false))
    end
    for i in padded_n+δn+start+1:2padded_n+1
        x = comm(basis, δn+start, i)==0x1
        z = comm(basis, padded_n+δn+start, i)==0x1
        z && mul_left!(basis, i, δn+start; phases=Val(false))
        x && mul_left!(basis, i, padded_n+δn+start; phases=Val(false))
    end
end

"""Perform the Symplectic Gram-Schmidt procedure that gives a Clifford operator canonically related to a given Pauli operator.

The algorithm is detailed in [koenig2014efficiently](@cite).

```jldoctest
julia> symplecticGS(P"X", padded_n=3)
X₁ ⟼ + X__
X₂ ⟼ + _X_
X₃ ⟼ + __X
Z₁ ⟼ + Z__
Z₂ ⟼ + _Z_
Z₃ ⟼ + __Z

julia> symplecticGS(P"Z")
X₁ ⟼ + Z
Z₁ ⟼ + X
```

See also: [`enumerate_cliffords`](@ref), [`clifford_cardinality`](@ref)."""
function symplecticGS(pauli::PauliOperator; padded_n=nqubits(pauli))
    n = nqubits(pauli)
    basis = zero(Tableau, 2padded_n+1, padded_n)
    δn = padded_n-n
    # fillup the padded tableau
    for i in 1:δn
        basis[i,i] = (true,false)
        basis[padded_n+i,i] = (false,true)
    end
    for i in δn+1:padded_n-1
        basis[i+1,i] = (true,false)
        basis[padded_n+i+1,i] = (false,true)
    end
    for i in 1:n basis[δn+1,δn+i] = pauli[i] end
    basis[padded_n+δn+1,padded_n] = (true,false)
    basis[2padded_n+1,padded_n] = (false,true)
    # end fillup
    doneupto = 1
    while doneupto <= n
        i = _findanticommGS(basis, doneupto, n, padded_n, δn)
        rowswap!(basis, padded_n+δn+doneupto, i; phases=Val(false))
        _eliminateGS(basis, doneupto, n, padded_n, δn)
        doneupto += 1
    end
    CliffordOperator((@view basis[1:2padded_n]))
end

@inline function int_to_bits(n,i)
    Bool[i>>d&0x1 for d in 0:n-1]
end

"""Give the i-th n-qubit Clifford operation, where i∈{1..2ⁿⁿΠⱼ₌₁ⁿ(4ʲ-1)}

The algorithm is detailed in [koenig2014efficiently](@cite).

See also: [`symplecticGS`](@ref), [`clifford_cardinality`](@ref)."""
function enumerate_cliffords(n,i;padded_n=n,onlycoset=false)
    enumerate_cliffords_slow(n,i;padded_n,onlycoset)
end

"""The O(n^4) implementation from [koenig2014efficiently](@cite) -- their algorithm seems wrong as ⟨w'₁|wₗ⟩=bₗ which is not always zero."""
function enumerate_cliffords_slow(n,i;padded_n=n,onlycoset=false) # TODO implement the faster n^3 algorithm
    @assert n<32 # Assuming 64 bit system
    s = 2^(2n)-1
    k = (i-1)%s+1
    pauli = PauliOperator(int_to_bits(2n,k))
    op = symplecticGS(pauli;padded_n)
    basis = tab(op)
    δn = padded_n-n
    idivs = (i-1)÷s
    for j in 0:n-1 # first n bits # w₁ ← w₁+vⱼ
        iszero(idivs>>j&0x1) || mul_left!(basis, padded_n+δn+1, δn+j+1; phases=Val(false))
    end
    for j in n:2n-2 # following n-1 bits # w₁ ← w₁+wⱼ
        iszero(idivs>>j&0x1) || mul_left!(basis, padded_n+δn+1, 2δn+j+2; phases=Val(false))
    end
    for j in 1:n-1 # first n-1 bits after first bit # wⱼ ← wⱼ+v₁ # missing from koenig2014efficiently
        iszero(idivs>>j&0x1) || mul_left!(basis, padded_n+δn+j+1, δn+1; phases=Val(false)) # δn+j+1+n
    end
    for j in n:2n-2 # following n-1 bits # vⱼ ← vⱼ+v₁ # missing from koenig2014efficiently
        iszero(idivs>>j&0x1) || mul_left!(basis, δn+j+2-n, δn+1; phases=Val(false))
    end
    if n==1 || onlycoset
        return op
    else
        subop = enumerate_cliffords(n-1,idivs÷2^(2n-1)+1;padded_n)
        return apply!(subop, op; phases=false)
    end
end

"""Give all n-qubit Clifford operations.

The algorithm is detailed in [koenig2014efficiently](@cite).

See also: [`symplecticGS`](@ref), [`clifford_cardinality`](@ref)."""
function enumerate_cliffords(n;padded_n=n,onlycoset=false) # TODO make an enumerator that does not compute the cosets from scratch each time
    if onlycoset
        cosetsize = (2^(2n)-1)*2^(2n-1)
        (enumerate_cliffords(n,i;padded_n,onlycoset=true) for i in 1:cosetsize)
    else
        (enumerate_cliffords(n,i;padded_n,onlycoset=false) for i in 1:clifford_cardinality(n;phases=false))
    end
end

@inline function _change_phases!(op, phases)
    tab(op).phases .= 0x2 .* phases
    op
end

"""Given an operator, return all operators that have the same tableau but different phases.

```jldoctest
julia> length(collect(enumerate_phases(tCNOT)))
16
```

See also: [`enumerate_cliffords`](@ref), [`clifford_cardinality`](@ref)."""
function enumerate_phases(op::CliffordOperator)
    n = nqubits(op)
    (_change_phases!(copy(op), int_to_bits(2n,i)) for i in 0:2^(2n)-1)
end

"""Given a set of operators, return all operators that have the same tableaux but different phases.

```jldoctest
julia> length(collect(enumerate_phases(enumerate_cliffords(2))))
11520
```

See also: [`enumerate_cliffords`](@ref), [`clifford_cardinality`](@ref)."""
function enumerate_phases(ops::Union{AbstractVector,Base.Generator})
    Iterators.flatten(Iterators.map(enumerate_phases, ops))
end
