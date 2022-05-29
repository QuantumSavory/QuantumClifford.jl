function halfclip!(state::AbstractStabilizer; phases::Bool=true)
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
            if k1!=columns
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
    end
    state
end