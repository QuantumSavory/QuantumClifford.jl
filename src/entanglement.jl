"""
Get the clipped gauge of a stablizer state. 
"""
function clip!(state::AbstractStabilizer; phases::Bool=true)
    _phases = Val(phases)
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
        k1 = findfirst(e->e&jsmall!=zerobit, xs[jbig,i:end] .| zs[jbig,i:end])
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
                        mul_left!(state, m, i; phases=_phases)
                    elseif (xs[jbig,m]⊻xs[jbig,i+1])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,i+1])&jsmall==zerobit
                        mul_left!(state, m, i+1; phases=_phases)
                    elseif (xs[jbig,m]⊻xs[jbig,i]⊻xs[jbig,i+1])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,i]⊻zs[jbig,i+1])&jsmall==zerobit
                        mul_left!(state, m, i; phases=_phases)
                        mul_left!(state, m, i+1; phases=_phases)
                    end
                end
                i += 2
            else # can only find k1
                # move k1 up to i
                rowswap!(state, k1, i; phases=phases)         
                # use it to eliminate others
                for m in i+1:rows
                    if (xs[jbig,m]⊻xs[jbig,i])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,i])&jsmall==zerobit
                        mul_left!(state, m, i; phases=_phases)
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
        k1 = findfirst(e->e&jsmall!=zerobit, xs[jbig,unfrozen_rows] .| zs[jbig,unfrozen_rows])
        
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
                        mul_left!(state, m, k1_row; phases=_phases)
                    end
                end
                # for other rows, use both
                for m in unfrozen_rows[k2+1:end]
                    if (xs[jbig,m]⊻xs[jbig,k1_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k1_row])&jsmall==zerobit
                        mul_left!(state, m, k1_row; phases=_phases)
                    elseif (xs[jbig,m]⊻xs[jbig,k2_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k2_row])&jsmall==zerobit
                        mul_left!(state, m, k2_row; phases=_phases)
                    elseif (xs[jbig,m]⊻xs[jbig,k1_row]⊻xs[jbig,k2_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k1_row]⊻zs[jbig,k2_row])&jsmall==zerobit
                        mul_left!(state, m, k1_row; phases=_phases)
                        mul_left!(state, m, k2_row; phases=_phases)
                    end
                end
                deleteat!(unfrozen_rows, (k1, k2))
            else # can only find k1
                # use it to eliminate others
                for m in unfrozen_rows[k1+1:end]
                    if (xs[jbig,m]⊻xs[jbig,k1_row])&jsmall==zerobit && (zs[jbig,m]⊻zs[jbig,k1_row])&jsmall==zerobit
                        mul_left!(state, m, k1_row; phases=_phases)
                    end
                end
                deleteat!(unfrozen_rows, k1)
            end
        end
    end
    state
end


"""
Get bigram which contains the location of endpoints.
"""
function get_bigram(state::AbstractStabilizer; do_clip::Bool=true)
    if do_clip
        clip!(state)
    end
    xzs = stabilizerview(state).xzs
    xs = @view xzs[1:end÷2,:]
    zs = @view xzs[end÷2+1:end,:]
    Tme = eltype(xzs)
    lowbit = Tme(0x1)
    zerobit = Tme(0x0)
    rows, columns = size(stabilizerview(state))
    xorzs = xs .| zs
    bigram = zeros(Int, rows, 2)
    for i in 1:rows
        bigram[i, 1] = findfirst(j->(
            jbig = _div(Tme,j-1)+1;
            jsmall = lowbit<<_mod(Tme,j-1);
            xorzs[jbig,i]&jsmall!=zerobit), 1:columns)
        bigram[i, 2] = findlast(j->(
            jbig = _div(Tme,j-1)+1;
            jsmall = lowbit<<_mod(Tme,j-1);
            xorzs[jbig,i]&jsmall!=zerobit), 1:columns)
    end
    bigram
end


"""
Get bipartite entanglement of a contiguous subsystem.
"""
function entanglement_cont(state::AbstractStabilizer, subsystem_ends; do_clip::Bool=true)
    bigram = get_bigram(state; do_clip=do_clip)
    leftend, rightend = subsystem_ends
    0.5 * count(r->(leftend<=r[1]<=rightend)⊻(leftend<=r[2]<=rightend), eachrow(bigram))
end
