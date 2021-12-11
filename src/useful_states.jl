"""Prepare one or more Bell pairs (with optional phases).

```jldoctest
julia> bell()
+ XX
+ ZZ

julia> bell(2)
+ XX__
+ ZZ__
+ __XX
+ __ZZ

julia> bell((true, false))
- XX
+ ZZ

julia> bell([true, false, true, true])
- XX__
+ ZZ__
- __XX
- __ZZ
```
"""
function bell end

function bell()
    copy(S"XX
           ZZ")
end

function bell(n::Int)
    tensor_pow(bell(), n)
end

function bell(phase::Tuple{Bool, Bool})
    s = bell()
    s.phases[1] = phase[1] ? 0x2 : 0x0
    s.phases[2] = phase[2] ? 0x2 : 0x0
    s
end

function bell(phases::AbstractVector{Tuple{Bool, Bool}})
    ⊗((bell(t) for t in phases)...)
end

function bell(phases::AbstractVector{Bool})
    s = bell(length(phases)÷2)
    for i in 1:length(s)
        phases[i]
        s.phases[i] = phases[i] ? 0x2 : 0x0
        s[i], s.phases[i]
    end
    s
end

"""
Prepare a GHZ state of n qubits.

```jldoctest
julia> ghz()
+ XXX
+ ZZ_
+ _ZZ

julia> ghz(2)
+ XX
+ ZZ

julia> ghz(4)
+ XXXX
+ ZZ__
+ _ZZ_
+ __ZZ
```
"""
function ghz end

function ghz()
    copy(S"XXX
           ZZI
           IZZ")
end

function ghz(n::Int)
    s = zero(Stabilizer, n)
    for i in 1:n
        s[1,i] = (true,false)
    end
    for i in 2:n
        s[i,i] = (false,true)
        s[i,i-1] = (false,true)
    end
    s
end

# TODO document these explicitly
# TODO cluster states, toric code, planar code, other codes from python libraries