xz2str(x,z) = join(toletter[e] for e in zip(x,z))

function xz2str_limited(x,z, limit=50)
    tupl = collect(zip(x,z))
    n = length(tupl)
    if (ismissing(limit) || limit >= n)
        return xz2str(x, z)
    end
    padding = limit÷2
    return join(toletter[tupl[i]] for i in 1:padding) * "⋯" * join(toletter[tupl[i]] for i in (n-padding):n)
end

_show(io::IO, p::PauliOperator, limit=50) = print(io, ["+ ","+i","- ","-i"][p.phase[]+1]*xz2str_limited(xbit(p),zbit(p), limit))

function Base.show(io::IO, p::PauliOperator)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        _show(io, p, 10)
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, p, max(2, sz[2]-7))
    else
        _show(io, p, missing)
    end
end

function _show(io::IO, t::Tableau, limit=50, limit_vertical=20)
    padding = max(1,limit_vertical÷2-3)
    n = size(t,1)
    range = 1:n
    if (!ismissing(limit_vertical) && limit_vertical < n+3)
        range = [1:padding; missing; (n-padding+2):n]
    end
    for i in range
        if ismissing(i)
            print(io," ⋮\n")
        else
            _show(io, t[i], max(1,limit-7))
            i!=n && println(io)
        end
    end
end

function Base.show(io::IO, t::Tableau)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        r,q = size(t)
        print(io, "Tableaux $r×$q")
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, t, sz[2], sz[1])
    else
        _show(io, t, missing, missing)
    end
end

function Base.show(io::IO, s::Stabilizer)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        r,q = size(s)
        print(io, "Stabilizer $r×$q")
    else
        show(io, tab(s))
    end
end

function Base.show(io::IO, s::MixedStabilizer)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        r,q = size(stabilizerview(s))
        print(io, "MixedStabilizer $r×$q")
    else
        show(io, stabilizerview(s))
    end
end

function Base.show(io::IO, d::Destabilizer)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        r,q = size(stabilizerview(d))
        print(io, "Destablizer $r×$q")
    elseif get(io, :limit, false)
        h,w = displaysize(io)
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(min(w-9,size(d.tab,2)-4),0))
        _show(io, destabilizerview(d).tab, w, h÷2)
        println(io, "\n𝒮𝓉𝒶𝒷" * "━"^max(min(w-7,size(d.tab,2)-2),0))
        _show(io, stabilizerview(d).tab, w, h÷2)
    else
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-4,0))
        _show(io, destabilizerview(d).tab, missing, missing)
        println(io, "\n𝒮𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-2,0))
        _show(io, stabilizerview(d).tab, missing, missing)
    end
end

function Base.show(io::IO, d::MixedDestabilizer)
    r = rank(d)
    q = nqubits(d)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        print(io, "MixedDestablizer $r×$q")
    elseif get(io, :limit, false)
        h,w = displaysize(io)
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(min(w-9,size(d.tab,2)-4),0))
        _show(io, destabilizerview(d).tab, w, h÷4)
        if r != q
            println(io)
            println(io, "𝒳ₗ" * "━"^max(min(w-5,size(d.tab,2)),0))
            _show(io, logicalxview(d).tab, w, h÷4)
        end
        println(io)
        println(io, "𝒮𝓉𝒶𝒷" * "━"^max(min(w-7,size(d.tab,2)-2),0))
        _show(io, stabilizerview(d).tab, w, h÷4)
        if r != q
            println(io)
            println(io, "𝒵ₗ" * "━"^max(min(w-5,size(d.tab,2)),0))
            _show(io, logicalzview(d).tab, w, h÷4)
        end
    else
        println(io, "𝒟ℯ𝓈𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-4,0))
        _show(io, destabilizerview(d).tab, missing, missing)
        if r != q
            println(io)
            println(io, "𝒳ₗ" * "━"^max(size(d.tab,2),0))
            _show(io, logicalxview(d).tab, missing, missing)
        end
        println(io)
        println(io, "𝒮𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-2,0))
        _show(io, stabilizerview(d).tab, missing, missing)
        if r != q
            println(io)
            println(io, "𝒵ₗ" * "━"^max(size(d.tab,2)),0)
            _show(io, logicalzview(d).tab, missing, missing)
        end
    end
end

function _show(io::IO, c::CliffordOperator, limit=50, limit_vertical=20)
    n = nqubits(c)
    nwidth = Int(ceil(log10(n+1)))
    _limit = limit-nwidth-10
    range = 1:n
    if (!ismissing(limit_vertical) && limit_vertical < n)
        padding = limit_vertical÷4
        range = [1:padding-1; missing; (n-padding+2):n]
    end
    for i in range
        if ismissing(i)
            print(" ⋮\n")
            continue
        end
        print(io, "X"*digits_substr(i,nwidth)*" ⟼ ")
        _show(io, c.tab[i], _limit)
        println(io)
    end
    for i in range
        if ismissing(i)
            print(" ⋮\n")
            continue
        end
        print(io, "Z"*digits_substr(i,nwidth)*" ⟼ ")
        _show(io, c.tab[i+n], _limit)
        i!=n && println(io)
    end
end

function Base.show(io::IO, c::CliffordOperator)
    if get(io, :compact, false) | haskey(io, :typeinfo)
        q = nqubits(c)
        print(io, "CliffordOperator on $q qubits")
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, c, sz[2], sz[1])
    else
        _show(io, c, missing, missing)
    end
end
