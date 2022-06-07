"""
$TYPEDSIGNATURES

Fix the clipped gauge of a stabilizer (in place).

Assumes the input is a valid stabilizer (all operators commute and have
real phases). It permits identity generators. But it is required that the number of generators should be the same as the number of qubits. 

```jldoctest
julia> ghz = S"XXXX
               ZZII
               IZZI
               IIZZ";

julia> canonicalize_clip!(ghz)
+ XXXX
+ ZZ__
+ _ZZ_
+ __ZZ
```

If `phases=false` is set, the canonicalization does not track the phases
in the tableau. (A speedup is expected but not benchmarked.)

```jldoctest
julia> s = S"-ZX
              XX"
- ZX
+ XX

julia> canonicalize_clip!(copy(s))
+iY_
+ XX

julia> canonicalize_clip!(copy(s), phases=false)
- Y_
+ XX
```

Introduced in [nahum2017quantum](@cite), with a more detailed explanation of the algorithm in Appendix A of [li2019measurement](@cite)
"""
function canonicalize_clip!(state::AbstractStabilizer; phases::Bool=true)
    _canonicalize_clip!(state; phases=Val(phases))
end
function  _canonicalize_clip!(state::AbstractStabilizer; phases::Val{B}=Val(true)) where B
    xzs = stabilizerview(state).xzs
    xs = @view xzs[1:end÷2,:]
    zs = @view xzs[end÷2+1:end,:]
    Tme = eltype(xzs)
    lowbit = Tme(0x1)
    zerobit = Tme(0x0)
    rows, columns = size(stabilizerview(state))
    # step 1: pregauge
    i = 1 # index to place used stab
    for j in 1:columns
        jbig = _div(Tme,j-1)+1
        jsmall = lowbit<<_mod(Tme,j-1)
        # find first row that is not I in col j
        k1 = findfirst(k->(xs[jbig,k] .| zs[jbig,k])&jsmall!=zerobit, i:rows)
        # find second row that is not I and not same as k1
        if k1!==nothing
            k1 += i-1
            k2 = findfirst(k->
                    jsmall & # take the bit
                    ((xs[jbig,k] | zs[jbig,k]) & # not identity
                    ((xs[jbig,k]⊻xs[jbig,k1]) | (zs[jbig,k]⊻zs[jbig,k1]))) != zerobit, # not same as k1
                k1+1:columns)
            if k2!==nothing
                k2 += k1
                # move k1 and k2 up to i and i+1
                rowswap!(state, k1, i; phases=phases)
                rowswap!(state, k2, i+1; phases=phases)
                # use them to eliminate others
                for m in i+2:rows
                    if (xs[jbig,m]⊻xs[jbig,i])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,i])&jsmall==zerobit
                        mul_left!(state, m, i; phases=phases)
                    elseif (xs[jbig,m]⊻xs[jbig,i+1])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,i+1])&jsmall==zerobit
                        mul_left!(state, m, i+1; phases=phases)
                    elseif (xs[jbig,m]⊻xs[jbig,i]⊻xs[jbig,i+1])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,i]⊻zs[jbig,i+1])&jsmall==zerobit
                        mul_left!(state, m, i; phases=phases)
                        mul_left!(state, m, i+1; phases=phases)
                    end
                end
                i += 2
            else # can only find k1
                # move k1 up to i
                rowswap!(state, k1, i; phases=phases)         
                # use it to eliminate others
                for m in i+1:rows
                    if (xs[jbig,m]⊻xs[jbig,i])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,i])&jsmall==zerobit
                        mul_left!(state, m, i; phases=phases)
                    end
                end
                i += 1
            end
        end
    end
    # step 2: gauge
    unfrozen_rows = Array(rows:-1:1)
    for j in columns:-1:1 # in reversed order to keep left ends
        jbig = _div(Tme,j-1)+1
        jsmall = lowbit<<_mod(Tme,j-1)
        # find first row that is not I in col j
        k1 = findfirst(k->(xs[jbig,k] .| zs[jbig,k])&jsmall!=zerobit, unfrozen_rows)
        
        # find second row that is not I and not same as k1
        if k1!==nothing
            k1_row = unfrozen_rows[k1]
            k2 = findfirst(k->
                    jsmall & # take the bit
                    ((xs[jbig,k] | zs[jbig,k]) & # not identity
                    ((xs[jbig,k]⊻xs[jbig,k1_row]) | (zs[jbig,k]⊻zs[jbig,k1_row]))) != zerobit, # not same as k1
                unfrozen_rows[k1+1:end])
            
            if k2!==nothing
                k2 += k1
                k2_row = unfrozen_rows[k2]
                # use them to eliminate others
                # for rows between k1 and k2, use k1 
                for m in unfrozen_rows[k1+1:k2-1]
                    if (xs[jbig,m]⊻xs[jbig,k1_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k1_row])&jsmall==zerobit
                        mul_left!(state, m, k1_row; phases=phases)
                    end
                end
                # for other rows, use both
                for m in unfrozen_rows[k2+1:end]
                    if (xs[jbig,m]⊻xs[jbig,k1_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k1_row])&jsmall==zerobit
                        mul_left!(state, m, k1_row; phases=phases)
                    elseif (xs[jbig,m]⊻xs[jbig,k2_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k2_row])&jsmall==zerobit
                        mul_left!(state, m, k2_row; phases=phases)
                    elseif (xs[jbig,m]⊻xs[jbig,k1_row]⊻xs[jbig,k2_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k1_row]⊻zs[jbig,k2_row])&jsmall==zerobit
                        mul_left!(state, m, k1_row; phases=phases)
                        mul_left!(state, m, k2_row; phases=phases)
                    end
                end
                deleteat!(unfrozen_rows, (k1, k2))
            else # can only find k1
                # use it to eliminate others
                for m in unfrozen_rows[k1+1:end]
                    if (xs[jbig,m]⊻xs[jbig,k1_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k1_row])&jsmall==zerobit
                        mul_left!(state, m, k1_row; phases=phases)
                    end
                end
                deleteat!(unfrozen_rows, k1)
            end
        end
    end
    state
end


"""
Get bigram which is a tableau containing the location of endpoints. Each row for a qubit, and each column for left (1) and right (2).

If `clip=false` is set, the function will calculate bigram without fixing clipped gauge. This is for the case where the input state is already in clipped gauge.

Introduced in [nahum2017quantum](@cite), with a more detailed explanation of the algorithm in Appendix A of [li2019measurement](@cite)

See also: [`canonicalize_clip!`](@ref)
"""
function bigram(state::AbstractStabilizer; clip::Bool=true)
    if clip
        clipped_state = canonicalize_clip!(copy(state))
    else
        clipped_state = state
    end
    xzs = stabilizerview(clipped_state).xzs
    xs = @view xzs[1:end÷2,:]
    zs = @view xzs[end÷2+1:end,:]
    Tme = eltype(xzs)
    lowbit = Tme(0x1)
    zerobit = Tme(0x0)
    rows, columns = size(stabilizerview(clipped_state))
    xorzs = xs .| zs
    bg = zeros(Int, rows, 2)
    for i in 1:rows
        bg[i, 1] = findfirst(j->(
            jbig = _div(Tme,j-1)+1;
            jsmall = lowbit<<_mod(Tme,j-1);
            xorzs[jbig,i]&jsmall!=zerobit), 1:columns)
        bg[i, 2] = findlast(j->(
            jbig = _div(Tme,j-1)+1;
            jsmall = lowbit<<_mod(Tme,j-1);
            xorzs[jbig,i]&jsmall!=zerobit), 1:columns)
        #TODO: check whether findfirst(j->|(tab[i,j]...), 1:columns) is faster
    end
    bg
end


"""
Get bipartite entanglement entropy of a subsystem, which is defined as entropy of the reduced density matrix.

Two backends are supported: `:clipping` uses clipping algorithm and supports only a contiguous subsystem, `:graph` uses an algorithm based on conversion to graph states.
"""
function entanglement_entropy end


"""
Get bipartite entanglement entropy of a contiguous subsystem by clipping algorithm.

If `clip=false` is set, the function will calculate entanglement entropy (by bigram) without fixing clipped gauge. This is for the case where the input state is already in clipped gauge.

See also: [`bigram`](@ref), [`canonicalize_clip!`](@ref)
"""
function entanglement_entropy(state::AbstractStabilizer, subsystem_range::UnitRange, ::Val{:clipping}; clip::Bool=true)
    bg = bigram(state; clip=clip)
    count(r->(r[1] in subsystem_range)⊻(r[2] in subsystem_range), eachrow(bg)) ÷ 2
end


import Graphs, Nemo, LinearAlgebra


"""
Get bipartite entanglement entropy by first converting the state to a graph.

Based on [hein2006entanglement](@cite).
"""
function entanglement_entropy(state::AbstractStabilizer, subsystem::AbstractVector, ::Val{:graph})
    graph = Graphs.Graph(state)
    adjmat = Graphs.adjacency_matrix(graph)
    other_subsystem = filter(i->!(i in collect(subsystem)), 1:Graphs.nv(graph))
    subadjmat = Nemo.matrix(Nemo.ResidueRing(Nemo.ZZ, 2), collect(adjmat[subsystem,other_subsystem]))
    LinearAlgebra.rank(subadjmat)
end


#TODO: function entanglement_entropy(state, range::AbstractVector, ::Val{:traceout})
