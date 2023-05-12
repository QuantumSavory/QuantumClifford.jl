_show(io::IO, p::PauliOperator, limit=50) = print(io, ["+ ","+i","- ","-i"][p.phase[]+1]*xz2str_limited(xbit(p),zbit(p), limit))

function Base.show(io::IO, p::PauliOperator)
    if get(io, :compact, false)
        _show(io, p, 10)
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, p, sz[2]-7)
    else
        _show(io, p, -1)
    end
end

function _show(io::IO, t::Tableau, limit=50, limit_vertical=20)
    padding = limit_vertical÷2-3
    n = size(t,1)
    range = 1:n
    if (limit_vertical < n && limit_vertical != -1)
        range = [1:padding; -1; (n-padding):n]
    end
    for i in range
        if (i == -1)
            print(io," ⋮\n")
            continue
        end
        _show(io, t[i], limit-7)
        i!=n && println(io)
    end
end

function Base.show(io::IO, t::Tableau)
    if get(io, :compact, false)
        r,q = size(t)
        print(io, "Tableaux $r×$q")
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, t, sz[2], sz[1])
    else
        _show(io, t, -1, -1)
    end
end

function Base.show(io::IO, s::Stabilizer)
    if get(io, :compact, false)
        r,q = size(s)
        print(io, "Stabilizer $r×$q")
    else
        show(io, tab(s))
    end
end

function Base.show(io::IO, s::MixedStabilizer)
    if get(io, :compact, false)
        r,q = size(stabilizerview(s))
        print(io, "MixedStabilizer $r×$q")
    else
        show(io, stabilizerview(s))
    end
end

function Base.show(io::IO, d::Destabilizer)
    if get(io, :compact, false)
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
        _show(io, destabilizerview(d).tab, -1, -1)
        println(io, "\n𝒮𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-2,0))
        _show(io, stabilizerview(d).tab, -1, -1)
    end
end

function Base.show(io::IO, d::MixedDestabilizer)
    r = rank(d)
    q = nqubits(d)
    if get(io, :compact, false)
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
        _show(io, destabilizerview(d).tab, -1, -1)
        if r != q
            println(io)
            println(io, "𝒳ₗ" * "━"^max(size(d.tab,2),0))
            _show(io, logicalxview(d).tab, -1, -1)
        end
        println(io)
        println(io, "𝒮𝓉𝒶𝒷" * "━"^max(size(d.tab,2)-2,0))
        _show(io, stabilizerview(d).tab, -1, -1)
        if r != q
            println(io)
            println(io, "𝒵ₗ" * "━"^max(size(d.tab,2)),0)
            _show(io, logicalzview(d).tab, -1, -1)
        end
    end
end

function _show(io::IO, c::CliffordOperator, limit=50, limit_vertical=20)
    n = nqubits(c)
    nwidth = Int(ceil(log10(n+1)))
    _limit = limit==-1 ? -1 : limit-nwidth-10
    range = 1:n
    if (limit_vertical < n && limit_vertical != -1)
        padding = limit_vertical÷4
        range = [1:padding-1; -1; (n-padding+2):n]
    end
    for i in range
        if (i == -1)
            print(" ⋮\n")
            continue
        end
        print(io, "X"*digits_substr(i,nwidth)*" ⟼ ")
        _show(io, c.tab[i], _limit)
        println(io)
    end
    for i in range
        if (i == -1)
            print(" ⋮\n")
            continue
        end
        print(io, "Z"*digits_substr(i,nwidth)*" ⟼ ")
        _show(io, c.tab[i+n], _limit)
        i!=n && println(io)
    end
end

function Base.show(io::IO, c::CliffordOperator)
    if get(io, :compact, false)
        q = nqubits(c)
        print(io, "Clifford $q qubits")
    elseif get(io, :limit, false)
        sz = displaysize(io)
        _show(io, c, sz[2], sz[1])
    else
        _show(io, c, -1, -1)
    end
end
