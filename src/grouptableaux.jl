using Graphs
using LinearAlgebra
"""
input: (generating) set of PauliOperator
output: Tableau
```jldoctest
julia> groupify(S"XZ XZ")
+ __
+ XZ
+ ZX
+ YY
```
"""
function groupify(s::Stabilizer) # add checks?
    n = length(s)
    group = zero(Stabilizer, 2^n, nqubits(s))
    for i in 0:2^n-1
        for (digit_order, j) in enumerate(digits(i, base=2, pad=n))
            if j == 1
                group[i+1] *= s[digit_order]
            end
        end
    end
    return group # TODO use Tableau instead of Stabilizer
end


"""
For a given non-necessarily-minimal generating set represented as a tableau `s`,
return the minimal generating set. Only real phases are permitted.
"""
function get_generating_set(s) # TODO potentially rename to `generator_tableau` or `generating_set`
    s, _, r = canonicalize!(copy(s), ranks=true)
    if r == 0
        gs = zero(Stabilizer, 1, nqubits(s))
        if (1im * zero(Stabilizer, 1, nqubits(s))[1]) in s
            gs[1] = 1im * gs[1]
        elseif (-1 * zero(Stabilizer, 1, nqubits(s))[1]) in s
            gs[1] = -1 * gs[1]
        end
        return gs
    else
        return s[1:r, :]
    end
end

"""
For a set of logical operators, rewrite in a set of logical x and logical z with the appropriate pairwise anticommutation
""" 
"""
Returns the full pauli group of the given length. Phases besides + are ignored by default, but can be included by setting phases = true. 
"""
function pauligroup(n, phases=false)
    if phases
        s = zero(Stabilizer, 4^(n + 1), n)
        paulis = ((false, false), (true, false), (false, true), (true, true))
        for (i, P) in enumerate(Iterators.product(Iterators.repeated(paulis, n)...))
            for (j, p) in enumerate(P)
                s[i, j] = p
            end
        end
        for i in 1:4^n
            s[i+4^n] = -1 * s[i]
        end
        for i in 4^n:2*4^n
            s[i+4^n] = -1im * s[i]
        end
        for i in 2*4^n:3*4^n
            s[i+4^n] = -1 * s[i]
        end
    end
    if !phases
        s = zero(Stabilizer, 4^n, n)
        paulis = ((false, false), (true, false), (false, true), (true, true))
        for (i, P) in enumerate(Iterators.product(Iterators.repeated(paulis, n)...))
            for (j, p) in enumerate(P)
                s[i, j] = p
            end
        end
    end
    return s
end

"""
For a given s, find all Pauli elements of the number of qubits that commute with all elements in s
"""
function normalize(s::Stabilizer)
    n = nqubits(s)
    pgroup = pauligroup(n)
    ptype = typeof(s[1])
    normalizer = ptype[]

    for p in pgroup
        commutes = 0
        for q in s
            if comm(p, q) == 0x01
                commutes = 1
            end
        end
        if commutes == 0
            push!(normalizer, p)
        end
    end

    return Stabilizer(normalizer)
end

"""
For a given s, returns the subset of elements in s that commute with all elements in s.
"""
function center(s::Stabilizer)  # hilariously inefficient. centerify? centrify?
    center = []
    for P in s
        commutes = 0
        for Q in s
            if comm(P, Q) == 0x01
                commutes = 1
                break
            end
        end
        if commutes == 0
            push!(center, P)
        end
    end
    stabilizer = zero(Stabilizer, length(center), nqubits(s))
    for i in eachindex(center)
        stabilizer[i] = center[i]
    end
    return stabilizer
end

function contract(s::PauliOperator, subset)
    if s[subset] == zero(PauliOperator, length(subset))  # more efficient way to do this?
        return s[setdiff(1:length(s), subset)]
    else
        return nothing
    end
end

function delete(s::PauliOperator, subset)
    return s[setdiff(1:length(s), subset)]
end

function delete(𝒮::Stabilizer, subset)
    return 𝒮[:, setdiff(1:nqubits(𝒮), subset)]
end



