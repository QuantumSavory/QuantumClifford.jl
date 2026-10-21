"""
$TYPEDEF

Minimum-weight Viterbi decoder for a [`TailBiting`](@ref) convolutional code
[ollivier2003description](@cite), [forney2007convolutional](@cite).

The checks starting on frame `t` act only on frames `t … t+ν`, so the decoder sweeps the frames
once, keeping for every possible error on the last `ν` frames (the trellis state, `4^(nν)` of them)
the lightest error consistent with the syndrome so far. The ring is closed by fixing the first `ν`
frames and requiring the sweep to wrap back onto them. The cost per syndrome is linear in the
number of frames `N`, against the exponential cost of a lookup table over all correctable errors.
Results are cached per syndrome.

```jldoctest convolutional
julia> tb = TailBiting(QuantumConvolutionalCode([[P"XXX", P"XZY"], [P"ZZZ", P"ZYX"]]), 6);

julia> d = ConvolutionalViterbiDecoder(tb);

julia> E = P"______X_______Z___"; s = comm(E, parity_checks(tb));

julia> g = QuantumClifford.ECC.decode(d, s); PauliOperator(0x0, g[1:18], g[19:end])
+ ______X_______Z___
```

See also: [`TailBiting`](@ref), [`TableDecoder`](@ref)
"""
struct ConvolutionalViterbiDecoder <: AbstractSyndromeDecoder
    """the code being decoded"""
    tb::TailBiting
    """its parity checks"""
    H::Stabilizer
    """number of Paulis on one frame, `4^n`"""
    nP::Int
    """weight of each frame Pauli"""
    weight::Vector{Int}
    """`bits[ℓ+1][e]`: mask of the generators whose part on offset `ℓ` anticommutes with frame Pauli `e`"""
    bits::Vector{Vector{UInt64}}
    """frame Paulis grouped by their mask on the last offset"""
    by_last::Dict{UInt64, Vector{Int}}
    """the frame Paulis, indexed like `weight`"""
    frame_paulis::Vector{PauliOperator}
    """decoded syndromes"""
    cache::Dict{Vector{Bool}, Vector{Bool}}
end

function ConvolutionalViterbiDecoder(tb::TailBiting)
    c = tb.code; n, ν, m = c.n, memory(c), length(c.generators)
    m <= 64 || throw(ArgumentError("`ConvolutionalViterbiDecoder` supports at most 64 generators per frame."))
    nP = 4^n
    ps = map(0:nP-1) do e
        p = zero(PauliOperator, n)
        for q in 1:n
            p[q] = (isodd(e >> (2q - 2)), isodd(e >> (2q - 1)))
        end
        p
    end
    weight = [count(q -> p[q] != (false, false), 1:n) for p in ps]
    bits = [[reduce(|, (UInt64(comm(c.generators[j][ℓ+1], p)) << (j - 1) for j in 1:m); init = UInt64(0)) for p in ps]
            for ℓ in 0:ν]
    by_last = Dict{UInt64, Vector{Int}}()
    for e in 1:nP
        push!(get!(by_last, bits[ν+1][e], Int[]), e)
    end
    ConvolutionalViterbiDecoder(tb, parity_checks(tb), nP, weight, bits, by_last, ps, Dict{Vector{Bool}, Vector{Bool}}())
end

parity_checks(d::ConvolutionalViterbiDecoder) = d.H

function decode(d::ConvolutionalViterbiDecoder, syndrome_sample)
    key = Vector{Bool}(syndrome_sample .!= 0)
    get!(() -> _viterbi_ring(d, key), d.cache, key)
end

# The trellis state is the tuple of Paulis on the last ν frames, stored as an integer in base 4^n
# (digit k = the frame k-1 places after the oldest one).
function _viterbi_ring(d::ConvolutionalViterbiDecoder, s::Vector{Bool})
    c, N = d.tb.code, d.tb.N
    n, ν, m, nP = c.n, memory(c), length(c.generators), d.nP
    masks = [reduce(|, (UInt64(s[(j-1)*N+t+1]) << (j - 1) for j in 1:m); init = UInt64(0)) for t in 0:N-1]
    nS = nP^ν
    INF = typemax(Int) ÷ 2
    best, best_frames = INF, Int[]
    cost, new = fill(INF, nS), fill(INF, nS)
    back = zeros(Int, nS, N)
    head = zeros(Int, ν)
    for start in 0:nS-1
        x = start
        for k in 1:ν
            head[k] = x % nP + 1; x ÷= nP
        end
        fill!(cost, INF)
        cost[start+1] = sum(d.weight[e] for e in head)
        for u in 0:N-1                 # choose frame u+ν, completing the checks that start on frame u
            fill!(new, INF)
            wrap = u + ν >= N
            for st in 0:nS-1
                cst = cost[st+1]
                cst >= INF && continue
                need, x = masks[u+1], st
                for k in 1:ν
                    need ⊻= d.bits[k][x % nP + 1]; x ÷= nP
                end
                cands = get(d.by_last, need, nothing)
                cands === nothing && continue
                for e in cands
                    wrap && e != head[u+ν-N+1] && continue
                    nst = (st ÷ nP) + (e - 1) * nP^(ν - 1)
                    nc = cst + (wrap ? 0 : d.weight[e])
                    if nc < new[nst+1]
                        new[nst+1] = nc; back[nst+1, u+1] = st
                    end
                end
            end
            cost, new = new, cost
        end
        if cost[start+1] < best
            best = cost[start+1]
            frames = zeros(Int, N + ν); st = start
            for u in N-1:-1:0
                frames[u+ν+1] = (st ÷ nP^(ν - 1)) % nP + 1
                st = back[st+1, u+1]
            end
            frames[1:ν] = head
            best_frames = frames[1:N]
        end
    end
    E = zero(PauliOperator, n * N)
    for t in 0:N-1, q in 1:n
        E[n*t+q] = d.frame_paulis[best_frames[t+1]][q]
    end
    stab_to_gf2(E)
end
