"""
$TYPEDEF

A quantum convolutional code, given by its *basic generators* repeated in every frame
[ollivier2003description](@cite), [forney2007convolutional](@cite).

Generator `j` acts on frames `t, t+1, …, t+ν` with the Paulis `generators[j][1], …, generators[j][ν+1]`,
each on the `n` qubits of one frame. Shorter generators are padded with identities up to the
memory `ν`. The code itself has no fixed length; use [`TailBiting`](@ref) to close it on a ring of
`N` frames and obtain an ordinary finite stabilizer code.

```jldoctest convolutional
julia> c = QuantumConvolutionalCode([[P"XXX", P"XZY"], [P"ZZZ", P"ZYX"]])
QuantumConvolutionalCode(n=3, m=2, ν=1: XXX + D·XZY, ZZZ + D·ZYX)

julia> rate(c)
1//3
```

### Fields
    $TYPEDFIELDS

See also: [`TailBiting`](@ref), [`ConvolutionalViterbiDecoder`](@ref)
"""
struct QuantumConvolutionalCode
    """qubits per frame"""
    n::Int
    """`generators[j][ℓ+1]` is the part of generator `j` on frame `t+ℓ`"""
    generators::Vector{Vector{PauliOperator}}
    function QuantumConvolutionalCode(gens::AbstractVector)
        isempty(gens) && throw(ArgumentError("A convolutional code needs at least one generator."))
        n = nqubits(first(first(gens)))
        all(nqubits(p) == n for g in gens for p in g) ||
            throw(ArgumentError("Every frame part of every generator must act on the same number of qubits."))
        ν = maximum(length.(gens)) - 1
        padded = [vcat(collect(g), [zero(PauliOperator, n) for _ in 1:ν+1-length(g)]) for g in gens]
        new(n, padded)
    end
end

"""Qubits per frame of a [`QuantumConvolutionalCode`](@ref)."""
frame_size(c::QuantumConvolutionalCode) = c.n

"""Memory `ν` of a [`QuantumConvolutionalCode`](@ref): a generator reaches `ν` frames beyond its first."""
memory(c::QuantumConvolutionalCode) = length(first(c.generators)) - 1

"""Encoded qubits per physical qubit in the stream, `(n - m)/n`, for `m` independent generators per frame."""
rate(c::QuantumConvolutionalCode) = (c.n - length(c.generators)) // c.n

"""The generators as Pauli polynomials in the delay operator `D`, e.g. `"XXX + D·XZY"`."""
function polynomial_form(c::QuantumConvolutionalCode)
    map(c.generators) do g
        terms = String[]
        for (ℓ, p) in enumerate(g)
            all(q -> p[q] == (false, false), 1:c.n) && continue
            s = join(("IXZY"[1 + p[q][1] + 2p[q][2]] for q in 1:c.n))
            push!(terms, ℓ == 1 ? s : ℓ == 2 ? "D·" * s : "D^$(ℓ-1)·" * s)
        end
        join(terms, " + ")
    end
end

"""
Whether every generator commutes with every shift of every generator, i.e. the shifted
symplectic product of the code vanishes [wilde2010convolutional](@cite). Only such codes define
a valid stabilizer code at every length.
"""
function shifted_commutes(c::QuantumConvolutionalCode)
    ν = memory(c)
    for g in c.generators, h in c.generators, s in -ν:ν
        acc = 0
        for ℓ in 0:ν
            0 <= ℓ + s <= ν || continue
            acc += comm(g[ℓ+1], h[ℓ+s+1])
        end
        isodd(acc) && return false
    end
    true
end

function Base.show(io::IO, c::QuantumConvolutionalCode)
    print(io, "QuantumConvolutionalCode(n=$(c.n), m=$(length(c.generators)), ν=$(memory(c)): ",
          join(polynomial_form(c), ", "), ")")
end

"""
$TYPEDEF

A [`QuantumConvolutionalCode`](@ref) closed into a ring of `N` frames (a tail-biting code
[forney2007convolutional](@cite)), as an ordinary stabilizer code on `n⋅N` qubits.

Qubit `q` of frame `t` (counting frames from 0) is qubit `n⋅t + q`, and check row `(j-1)⋅N + t + 1` is
generator `j` shifted to start on frame `t`.

```jldoctest convolutional
julia> tb = TailBiting(QuantumConvolutionalCode([[P"XXX", P"XZY"], [P"ZZZ", P"ZYX"]]), 4);

julia> code_n(tb), code_k(tb)
(12, 4)

julia> parity_checks(tb)
+ XXXXZY______
+ ___XXXXZY___
+ ______XXXXZY
+ XZY______XXX
+ ZZZZYX______
+ ___ZZZZYX___
+ ______ZZZZYX
+ ZYX______ZZZ
```

### Fields
    $TYPEDFIELDS

See also: [`ConvolutionalViterbiDecoder`](@ref)
"""
struct TailBiting <: AbstractQECC
    """the convolutional code"""
    code::QuantumConvolutionalCode
    """number of frames in the ring"""
    N::Int
    function TailBiting(code::QuantumConvolutionalCode, N::Int)
        N > memory(code) || throw(ArgumentError("The ring needs more than ν = $(memory(code)) frames."))
        shifted_commutes(code) || throw(ArgumentError("The shifted generators of the code do not commute."))
        new(code, N)
    end
end

iscss(::Type{TailBiting}) = nothing

code_n(tb::TailBiting) = tb.code.n * tb.N

function parity_checks(tb::TailBiting)
    c, N, n = tb.code, tb.N, tb.code.n
    rows = PauliOperator[]
    for g in c.generators, t in 0:N-1
        p = zero(PauliOperator, n * N)
        for (ℓ, part) in enumerate(g), q in 1:n
            part[q] == (false, false) && continue
            i = n * mod(t + ℓ - 1, N) + q
            p[i] = (p[i][1] ⊻ part[q][1], p[i][2] ⊻ part[q][2])
        end
        push!(rows, p)
    end
    Stabilizer(rows)
end
