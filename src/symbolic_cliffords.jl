using Random: AbstractRNG, GLOBAL_RNG

abstract type AbstractSymbolicOperator <: AbstractCliffordOperator end
abstract type AbstractSingleQubitOperator <: AbstractSymbolicOperator end
abstract type AbstractTwoQubitOperator <: AbstractSymbolicOperator end

"""A "symbolic" Hadamard operator which permits faster multiplication than an operator expressed as an explicit tableau.

```jldoctest
julia> sh = sHadamard(1)
Symbolic single-qubit gate on qubit 1
X ⟼ + Z
Z ⟼ + X

julia> typeof(sh)
sHadamard

julia> h = CliffordOperator(sh) # Transforming it back into an explicit tableau representation
X ⟼ + Z
Z ⟼ + X

julia> typeof(h)
CliffordOperator{Vector{UInt8}, LinearAlgebra.Adjoint{UInt64, Matrix{UInt64}}}
```

Notice that in tensor products and other similar operations, this symbolic operator gets
transformed into a single-qubit tableau, losing the `qubit` index information.

See also: [`SingleQubitOperator`](@ref)
"""
struct sHadamard <: AbstractSingleQubitOperator
    q::Int
end

"""A "symbolic" Phase operator which permits faster multiplication than an operator expressed as an explicit tableau.

```jldoctest
julia> sp = sPhase(1)
Symbolic single-qubit gate on qubit 1
X ⟼ + Y
Z ⟼ + Z

julia> typeof(sp)
sPhase

julia> p = CliffordOperator(sp) # Transforming it back into an explicit tableau representation
X ⟼ + Y
Z ⟼ + Z

julia> typeof(p)
CliffordOperator{Vector{UInt8}, LinearAlgebra.Adjoint{UInt64, Matrix{UInt64}}}
```

Notice that in tensor products and other similar operations, this symbolic operator gets
transformed into a single-qubit tableau, losing the `qubit` index information.

See also: [`SingleQubitOperator`](@ref)
"""
struct sPhase <: AbstractSingleQubitOperator
    q::Int
end

"""A "symbolic" single-qubit operator which permits faster multiplication than an operator expressed as an explicit tableau.

```jldoctest
julia> op = SingleQubitOperator(1, true, true, true, false, true, true) # Tableau components and phases
Symbolic single-qubit gate on qubit 1
X ⟼ - Y
Z ⟼ - X

julia> typeof(op)
SingleQubitOperator

julia> t_op = CliffordOperator(op) # Transforming it back into an explicit tableau representation
X ⟼ - Y
Z ⟼ - X

julia> typeof(t_op)
CliffordOperator{Vector{UInt8}, LinearAlgebra.Adjoint{UInt64, Matrix{UInt64}}}

julia> typeof(op ⊗ op) # Tensor and other operations are possible but some of them revert to tableau representation
CliffordOperator{Vector{UInt8}, LinearAlgebra.Adjoint{UInt64, Matrix{UInt64}}}
```

Notice that in tensor products and other similar operations, this symbolic operator gets
transformed into a single-qubit tableau, losing the `qubit` index information.

See also: [`sHadamard`](@ref), [`sPhase`](@ref), [`CliffordOperator`](@ref)
"""
struct SingleQubitOperator <: AbstractSingleQubitOperator
    q::Int
    xx::Bool
    xz::Bool
    zx::Bool
    zz::Bool
    px::Bool
    pz::Bool
end

SingleQubitOperator(h::sHadamard) = SingleQubitOperator(h.q, false, true, true, false, false, false)
SingleQubitOperator(p::sPhase) = SingleQubitOperator(p.q, true, true, false, true, false, false)
function SingleQubitOperator(op::CliffordOperator, qubit)
    @assert nqubits(op)==1 "You are trying to convert a multiqubit `CliffordOperator` into a symbolic `SingleQubitOperator`."
    SingleQubitOperator(qubit,op.tab[1,1]...,op.tab[2,1]...,(~).(iszero.(op.tab.phases))...)
end
SingleQubitOperator(op::CliffordOperator) = SingleQubitOperator(op, 1)
CliffordOperator(op::AbstractSingleQubitOperator) = CliffordOperator(SingleQubitOperator(op))
CliffordOperator(op::SingleQubitOperator) = CliffordOperator(Stabilizer([op.px ? 0x2 : 0x0, op.pz ? 0x2 : 0x0],[op.xx op.xz; op.zx op.zz]))
function Base.show(io::IO, op::AbstractSingleQubitOperator)
    print(io, "Symbolic single-qubit gate on qubit $(op.q)\n")
    show(io, CliffordOperator(op))
end

@inline getshift(Tme::Type,col::Int) = _mod(Tme,col-1)
@inline getmask(Tme::Type,col::Int) = Tme(0x1)<<getshift(Tme,col)
@inline getbigindex(Tme::Type,col::Int) = _div(Tme,col-1)+1

Base.@propagate_inbounds function getxbit(s, r, c)
    Tme = eltype(s.xzs)
    s.xzs[r,getbigindex(Tme,c)]&getmask(Tme,c)
end
Base.@propagate_inbounds function getzbit(s, r, c)
    Tme = eltype(s.xzs)
    s.xzs[r,end÷2+getbigindex(Tme,c)]&getmask(Tme,c)
end
Base.@propagate_inbounds function setxbit(s, r, c, x)
    Tme = eltype(s.xzs)
    cbig = getbigindex(Tme,c)
    s.xzs[r,cbig] &= ~getmask(Tme,c)
    s.xzs[r,cbig] |= x
end
Base.@propagate_inbounds function setzbit(s, r, c, z)
    Tme = eltype(s.xzs)
    cbig = getbigindex(Tme,c)
    s.xzs[r,end÷2+cbig] &= ~getmask(Tme,c)
    s.xzs[r,end÷2+cbig] |= z
end
Base.@propagate_inbounds setxbit(s, r, c, x, shift) = setxbit(s, r, c, x<<shift)
Base.@propagate_inbounds setzbit(s, r, c, z, shift) = setzbit(s, r, c, z<<shift)

function apply!(stab::AbstractStabilizer, h::sHadamard; phases::Bool=true)
    s = tab(stab)
    c = h.q
    @inbounds @simd for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        setxbit(s, r, c, z)
        setzbit(s, r, c, x)
        phases && x!=0 && z!=0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function apply!(stab::AbstractStabilizer, p::sPhase; phases::Bool=true)
    s = tab(stab)
    c = p.q
    @inbounds @simd for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        #setxbit no op
        setzbit(s, r, c, x⊻z)
        phases && x!=0 && z!=0 && (s.phases[r] = (s.phases[r]+0x2)&3)
    end
    stab
end

function apply!(stab::AbstractStabilizer, op::SingleQubitOperator; phases::Bool=true) # TODO generated functions that simplify the whole `if phases` branch might be a good optimization
    s = tab(stab)
    c = op.q
    Tme = eltype(s.xzs)
    sh = getshift(Tme, c)
    xx,zx,xz,zz = Tme.((op.xx,op.zx,op.xz,op.zz)) .<< sh
    anticom = ~iszero((~zz & xz & ~xx & zx) | ( zz & ~xz & xx & zx) | (zz &  xz & xx & ~zx))
    @inbounds @simd for r in eachindex(s)
        x = getxbit(s, r, c)
        z = getzbit(s, r, c)
        setxbit(s, r, c, (x&xx)⊻(z&zx))
        setzbit(s, r, c, (x&xz)⊻(z&zz))
        if phases
            if op.px && ~iszero(x)
                s.phases[r] += 0x2
                s.phases[r] &= 3
            end
            if op.pz && ~iszero(z)
                s.phases[r] += 0x2
                s.phases[r] &= 3
            end
            if ~iszero(x&z) && anticom
                s.phases[r] += 0x2
                s.phases[r] &= 3
            end
        end
    end
    stab
end

const all_single_qubit_patterns = (
    (true, false, false, true),
    (false, true, true, true),
    (true, true, true, false),
    (false, true, true, false),
    (true, false, true, true),
    (true, true, false, true)
)

"""Generate a symbolic single-qubit gate given its index. Optionally, set non-trivial phases."""
function enumerate_single_qubit_gates(index; qubit=1, phases=nothing)
    @assert index<=6 "Only 6 single-qubit gates exit, up to the choice of phases"
    if isnothing(phases)
        if index==4
            return sHadamard(qubit)
        elseif index==6
            return sPhase(qubit)
        else
            return SingleQubitOperator(qubit, all_single_qubit_patterns[index]..., false, false)
        end
    else
        return SingleQubitOperator(qubit, all_single_qubit_patterns[index]..., phases...)
    end
end

"""Random symbolic single-qubit Clifford applied to qubit at index `qubit`.

See also: [`SingleQubitOperator`](@ref), [`random_clifford`](@ref)
"""
function random_clifford1(rng::AbstractRNG, qubit)
    return enumerate_single_qubit_gates(rand(rng,1:6),qubit=qubit,phases=rand(rng,Bool,2))
end
random_clifford1(qubit) = random_clifford1(GLOBAL_RNG, qubit)


"""A "symbolic" CNOT
"""
struct sCNOT <: AbstractSingleQubitOperator
    q1::Int
    q2::Int
end

"""A "symbolic" SWAP
"""
struct sSWAP <: AbstractSingleQubitOperator
    q1::Int
    q2::Int
end

function apply!(stab::AbstractStabilizer, cnot::sCNOT; phases::Bool=true)
    s = tab(stab)
    q1 = cnot.q1
    q2 = cnot.q2
    Tme = eltype(s.xzs)
    shift = getshift(Tme, q1) - getshift(Tme, q2)
    @inbounds @simd for r in eachindex(s)
        x1 = getxbit(s, r, q1)
        z1 = getzbit(s, r, q1)
        x2 = getxbit(s, r, q2)
        z2 = getzbit(s, r, q2)
        setxbit(s, r, q2, x2⊻(x1<<-shift))
        setzbit(s, r, q1, z1⊻(z2<<shift))
    end
    stab
end

function apply!(stab::AbstractStabilizer, swap::sSWAP; phases::Bool=true)
    s = tab(stab)
    q1 = swap.q1
    q2 = swap.q2
    Tme = eltype(s.xzs)
    shift = getshift(Tme, q1) - getshift(Tme, q2)
    @inbounds @simd for r in eachindex(s)
        x1 = getxbit(s, r, q1)
        z1 = getzbit(s, r, q1)
        x2 = getxbit(s, r, q2)
        z2 = getzbit(s, r, q2)
        setxbit(s, r, q1, x2, shift)
        setzbit(s, r, q1, z2, shift)
        setxbit(s, r, q2, x1, -shift)
        setzbit(s, r, q2, z1, -shift)
    end
    stab
end
