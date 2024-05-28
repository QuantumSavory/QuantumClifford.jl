#using QuantumClifford
#= using GLMakie
using GraphMakie =#
using Graphs
using LinearAlgebra
#using IterTools


#

# center
# normalizer
# `invert' arbitrary groups

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
return the minimal generating set.
"""
function get_generating_set(s) # TODO potentially rename to `generator_tableau` or `generating_set`
    s, _, r = canonicalize!(copy(s), ranks=true)
    return s[1:r,:]
end

"""
For a set of logical operators, rewrite in a set of logical x and logical z with the appropriate pairwise anticommutation
""" # TODO unfinished, including docs
function logical_operator_canonicalize(s)
    # TODO potentially use
    # julia> MixedDestabilizer(S"XXX") |> logicalxview
    # julia> MixedDestabilizer(S"XXX") |> logicalzview
    # or stabilizerview or destabilizerview
    # @assert length(s) == length(get_generating_set(groupify(s)))
    mx = MixedDestabilizer(s)
    return vcat(logicalxview(mx), logicalzview(mx), s), rank(mx) # TODO use Tableau instead of Stabilizer
end


##
# s = S"XZ ZX IY"
# s = S"ZZ XZ YY"
# println("")
# t = logical_operator_canonicalize(s)
# ##

function commutavise(s::Stabilizer) # TODO take a Tableau
    @assert length(s) == length(groupify(s))
    embedding = []
    # currently this works only on generating sets for abelian groups
    # maybe split in commutavise
    # news = ... with n+(n-rank) columns
    # news[:, 1:n] = s
    # add some value to news[:, n+1:2n-rank]
    for P in s
        for Q in s
            if comm(P, Q) == 0x01
                # finish
            end
        end
    end
end

# s = S"XZ ZX IY"
# s = groupify(s)
# println("")
# print(get_generating_set(s))
# ##

function pauligroup(n, phases = false) # ignores phases #TODO add arg to control whether phases are ignored
    s = zero(Stabilizer, 4^(n+1), n)
    paulis = ((false, false), (true, false), (false, true), (true, true))
    for (i,P) in enumerate(Iterators.product(Iterators.repeated(paulis, n)...))
        for (j,p) in enumerate(P)
            s[i,j] = p
        end
    end
    if phases
        for i in 1:4^n
            s[i+4^n] = -1* s[i]
        end
    end
    return s
end

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
    return Stabilizer(center)
end

# s = S"XXII IIXX ZIZI IZIZ"
# s = groupify(s)

#s = S"II YY"
#s = groupify(s)

#println(normalize(s))

##
# Pgroup = pauligroup(2)
# println(Pgroup)
##

# scenter = center(s)
# print(scenter)
# println(canonicalize!(s)[1:,])

# function center(s::Stabilizer) # centerify? centrify?
#     #  # returns pairs of non-commuting generators
# end

# function



function contract(s::PauliOperator, subset)
    if s[subset] == zero(PauliOperator, length(subset))  # more efficient way to do this?
        return s[setdiff(1:length(s), subset)]
    else
        return nothing
    end
end

function contract(𝒮::Stabilizer, subset)  # change mathcal notation since not a group?
    # TODO: contraction on a generating set does not yield a generating set of the contraction of <generating set>!
    return Stabilizer([P for P in contract.(𝒮, (subset,)) if !isnothing(P)])  # make more efficient
end

function delete(s::PauliOperator, subset)
    return s[setdiff(1:length(s), subset)]
end

function delete(𝒮::Stabilizer, subset)
    return 𝒮[:, setdiff(1:nqubits(𝒮), subset)]
end



# s1 = S"XIXIZIII ZZIIXIII ZIZIIXII IZIZIIXI ZIIIXIIX IIXIIZII IIIXIIZI IXIXZIIZ"
# s1 = groupify(s1)

# println("______!!!")
# s1 = delete(s1, [6, 7, 8])
# s1 = contract(s1, [5])


# s2 = S"XIXIZIII ZZIIXIII ZIZIIXII IZIZIIXI ZIIIXIIX IIXIIZII IIIXIIZI IXIXZIIZ"
# s2 = groupify(s2)

# println("______!!!")
# s2 = contract(s2, [6, 7, 8])
# s2 = contract(s2, [5])


# println(canonicalize!(s))
# println("______!!!")

# logical operators
# s = contract(s, [6, 7, 8])
# s = delete(s, [5])
#println(canonicalize!(s))


# s_contract = contract(s, [1,2])


# println("-----")

# s = S"XZ ZX"
# s = groupify(s)

# s = S"XZZ ZXZ ZZX"
# s = groupify(s)

# s = Stabilizer(wheel_graph(6))
# s = groupify(s)

# deleted = delete(s, [1])
# println(canonicalize!(deleted))
# println("_---")
# contracted = contract(s, [1])
# println(canonicalize!(contracted)[1:4])

# println("---")
# using QuantumClifford.ECC

# println(canonicalize!(parity_checks(Perfect5())[:,[1, 2, 3, 5, 4]]))
# println(typeof(groupify(S"XZ ZX")))


# # md = MixedDestabilizer(s); groupify(copy(vcat(stabilizerview(md), logicalxview(md), logicalzview(md)))) x
