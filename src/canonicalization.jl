"""
$TYPEDSIGNATURES

Canonicalize a stabilizer (in place).

Assumes the input is a valid stabilizer (all operators commute and have
real phases). It permits redundant generators and identity generators.

```jldoctest
julia> ghz = S"XXXX
               ZZII
               IZZI
               IIZZ";


julia> canonicalize!(ghz)
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ

julia> canonicalize!(S"XXXX
                       IZZI
                       IIZZ")
+ XXXX
+ _Z_Z
+ __ZZ
```

Not all rows in the tableau in the next example are independent:

```jldoctest
julia> canonicalize!(S"XXXX
                       ZZII
                       IZZI
                       IZIZ
                       IIZZ")
+ XXXX
+ Z__Z
+ _Z_Z
+ __ZZ
+ ____
```

In cases of lower rank, more advanced tableau structures might be better.
For instance the [`MixedStabilizer`](@ref) or [`MixedDestabilizer`](@ref)
structures (you can read more about them in the [Data Structures section](@ref Choosing-Appropriate-Data-Structure)
of the documentation).

If `phases=false` is set, the canonicalization does not track the phases
in the tableau, leading to significant (constant factor) speedup.

```jldoctest
julia> s = S"-ZX
              XZ"
- ZX
+ XZ

julia> canonicalize!(copy(s), phases=false)
- XZ
+ ZX

julia> canonicalize!(copy(s))
+ XZ
- ZX
```

If `ranks=true` is set, the last pivot indices for the X and Z stage of
the canonicalization are returned as well.

```jldoctest
julia> s = S"XXXX
             ZZII
             IZIZ
             ZIIZ";


julia> _, ix, iz = canonicalize!(s, ranks=true); ix, iz
(1, 3)

julia> s
+ XXXX
+ Z__Z
+ _Z_Z
+ ____
```

Based on [garcia2012efficient](@cite).

See also: [`canonicalize_rref!`](@ref), [`canonicalize_gott!`](@ref)
"""
function canonicalize!(state::AbstractStabilizer; phases::Bool=true, ranks::Bool=false)
    @valbooldispatch _canonicalize!(state; phases=Val(phases), ranks=Val(ranks)) phases ranks
end
function _canonicalize!(state::AbstractStabilizer; phases::Val{Bphases}=Val(true), ranks::Val{Branks}=Val(false)) where {Bphases,Branks}
    tab = stabilizerview(state)
    rows, columns = size(stabilizerview(state))
    i = 1
    for j in 1:columns
        # find first row with X or Y in col `j`
        k = findfirst(ii->tab[ii,j][1],i:rows)
        if k !== nothing
            k += i-1
            rowswap!(state, k, i; phases)
            for m in 1:rows
                if tab[m,j][1] && m!=i # if X or Y
                    mul_left!(state, m, i; phases)
                end
            end
            i += 1
        end
    end
    rx = i
    for j in 1:columns
        # find first row with Z in col `j`
        k = findfirst(ii->tab[ii,j][2],i:rows)
        if k !== nothing
            k += i-1
            rowswap!(state, k, i; phases)
            for m in 1:rows
                if tab[m,j][2] && m!=i # if Z or Y
                    mul_left!(state, m, i; phases)
                end
            end
            i += 1
        end
    end
    if Branks
        return state, rx-1, i-1
    else
        return state
    end
end

"""
$TYPEDSIGNATURES

Canonicalize a stabilizer (in place) along only some columns.

This uses different canonical form from [`canonicalize!`](@ref). It also indexes in
reverse in order to make its use in [`traceout!`](@ref) more efficient.
Its use in `traceout!` is its main application.

It returns the (in place) modified state and the index of the last pivot.

Based on [audenaert2005entanglement](@cite).

See also: [`canonicalize!`](@ref), [`canonicalize_gott!`](@ref)
"""
function canonicalize_rref!(state::AbstractStabilizer, colindices; phases::Bool=true)
    @valbooldispatch _canonicalize_rref!(state, colindices; phases=Val(phases)) phases
end
function _canonicalize_rref!(state::AbstractStabilizer, colindices; phases::Val{B}=Val(true)) where B
    tab = stabilizerview(state)
    rows, columns = size(stabilizerview(state))
    i = rows
    for j in colindices
        k = findfirst(ii->tab[ii,j][1],1:i)
        if k !== nothing
            rowswap!(state, k, i; phases)
            for m in 1:rows
                if tab[m,j][1] && m!=i # if X or Y
                    mul_left!(state, m, i; phases)
                end
            end
            i -= 1
        end
        k = findfirst(ii->tab[ii,j][2],1:i)
        if k !== nothing
            rowswap!(state, k, i; phases)
            for m in 1:rows
                if tab[m,j][2] && m!=i # if Z or Y
                    mul_left!(state, m, i; phases)
                end
            end
            i -= 1
        end
    end
    state, i
end

"""
$TYPEDSIGNATURES
"""
canonicalize_rref!(state::AbstractStabilizer; phases::Bool=true) = canonicalize_rref!(state, 1:nqubits(state); phases=phases)

function gott_standard_form_indices(chunks2D, rows, cols; skip=0)::Tuple{Vector{Int},Int}
    goodindices = Int[]
    j = 1
    r = 1
    for r in skip+1:rows
        i = unsafe_bitfindnext_(chunks2D[:,r],skip+1)
        isnothing(i) && break
        (i ∈ goodindices)::Bool && continue
        push!(goodindices, i)
    end
    rank = length(goodindices)
    if rank>0
        badindices = [r for r in 1+skip:goodindices[end] if !(r ∈ goodindices)]
        return vcat(1:skip, goodindices, badindices, goodindices[end]+1:cols), rank
    else
        return collect(1:cols), rank # without the collect it is not type stable; TODO is the collect making this slow?
    end
end

"""
Inplace Gottesman canonicalization of a tableau.

This uses different canonical form from [`canonicalize!`](@ref).
It is used in the computation of the logical X and Z operators
of a [`MixedDestabilizer`](@ref).

It returns the (in place) modified state, the indices of the last pivot
of both Gaussian elimination steps, and the permutations that have been used
to put the X and Z tableaux in standard form.

Based on [gottesman1997stabilizer](@cite).

See also: [`canonicalize!`](@ref), [`canonicalize_rref!`](@ref)
"""
function canonicalize_gott!(stabilizer::Stabilizer; phases::Bool=true)
    @valbooldispatch _canonicalize_gott!(stabilizer; phases=Val(phases)) phases
end
function _canonicalize_gott!(stabilizer::Stabilizer; phases::Val{B}=Val(true)) where {B}
    xzs = tab(stabilizer).xzs
    rows, columns = size(stabilizer)
    i = 1
    for j in 1:columns
        # find first row with X or Y in col `j`
        k = findfirst(ii->stabilizer[ii,j][1],i:rows)
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases)
            for m in 1:rows
                if stabilizer[m,j][1] && m!=i # if X or Y
                    mul_left!(stabilizer, m, i; phases)
                end
            end
            i += 1
        end
    end
    xperm, r = gott_standard_form_indices((@view xzs[1:end÷2,:]),rows,columns)
    permute!(stabilizer,xperm)
    i = r+1
    for j in r+1:columns
        # find first row with Z in col `j`
        k = findfirst(ii->stabilizer[ii,j][2],i:rows)
        if k !== nothing
            k += i-1
            rowswap!(stabilizer, k, i; phases)
            for m in 1:rows
                if stabilizer[m,j][2] && m!=i # if Z or Y
                    mul_left!(stabilizer, m, i; phases)
                end
            end
            i += 1
        end
    end
    zperm, s = gott_standard_form_indices((@view xzs[end÷2+1:end,:]),rows,columns,skip=r)
    permute!(stabilizer,zperm)
    stabilizer, r, s, xperm, zperm
end
