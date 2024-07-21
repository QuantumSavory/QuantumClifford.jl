"""
Generate a Pauli operator by using operators from a given the Stabilizer.

**It assumes the stabilizer is already canonicalized.** It modifies
the Pauli operator in place, generating it in reverse, up to a phase.
That phase is left in the modified operator, which should be the identity up to a phase.
Returns the new operator and the list of indices denoting the elements of
`stabilizer` that were used for the generation.

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

julia> generate!(P"-ZIZI", ghz)
(- ____, [2, 4])
```

When the Pauli operator can not be generated by the given tableau, `nothing` is returned.

```jldoctest
julia> generate!(P"XII",canonicalize!(S"ZII")) === nothing
true

julia> generate!(P"XII",canonicalize!(S"XII")) === nothing
false
```
"""
function generate!(pauli::PauliOperator, stabilizer::Stabilizer; phases::Bool=true, saveindices::Bool=true)
    @valbooldispatch _generate!(pauli, stabilizer; phases=Val(phases), saveindices=Val(saveindices)) phases saveindices
end

function _generate!(pauli::PauliOperator{Tz,Tv}, stabilizer::Stabilizer{Tableau{Tzv,Tm}}; phases::Val{PHASES}=Val(true), saveindices::Val{SAVEIDX}=Val(true)) where {Tz<:AbstractArray{UInt8,0}, Tzv<:AbstractVector{UInt8}, Tme<:Unsigned, Tv<:AbstractVector{Tme}, Tm<:AbstractMatrix{Tme}, PHASES, SAVEIDX} # TODO there is stuff that can be abstracted away here and in canonicalize!
    xzs = tab(stabilizer).xzs
    xs = @view xzs[1:end÷2,:]
    zs = @view xzs[end÷2+1:end,:]
    lowbit = Tme(0x1)
    zerobit = Tme(0x0)
    px,pz = xview(pauli), zview(pauli)
    used_indices = Int[]
    used = 0
    # remove Xs
    while (i=unsafe_bitfindnext_(px,1); i !== nothing) # TODO awkward notation due to https://github.com/JuliaLang/julia/issues/45499
        jbig = _div(Tme,i-1)+1
        jsmall = lowbit<<_mod(Tme,i-1)
        candidate = findfirst(e->e&jsmall!=zerobit, # TODO some form of reinterpret might be faster than equality check
                              xs[jbig,used+1:end])
        if isnothing(candidate)
            return nothing
        else
            used += candidate
        end
        mul_left!(pauli, stabilizer, used, phases=phases)
        SAVEIDX && push!(used_indices, used)
    end
    # remove Zs
    while (i=unsafe_bitfindnext_(pz,1); i !== nothing) # TODO awkward notation due to https://github.com/JuliaLang/julia/issues/45499
        jbig = _div(Tme,i-1)+1
        jsmall = lowbit<<_mod(Tme,i-1)
        candidate = findfirst(e->e&jsmall!=zerobit, # TODO some form of reinterpret might be faster than equality check
                              zs[jbig,used+1:end])
        if isnothing(candidate)
            return nothing
        else
            used += candidate
        end
        mul_left!(pauli, stabilizer, used, phases=phases)
        SAVEIDX && push!(used_indices, used)
    end
    all(iszero, pauli.xz) || return nothing # Need to check due to cases like generate!(P"_Z", S"XZ") that bypass the X checks # TODO do this better, without this extra loop through p.xz
    if SAVEIDX
        return pauli, used_indices
    else
        return pauli
    end
end

"""
$TYPEDSIGNATURES

Project the state of a Stabilizer on the two eigenspaces of a Pauli operator.

Assumes the input is a valid stabilizer. The projection is done inplace on
that stabilizer and it does not modify the projection operator.

It returns

 - a stabilizer that might not be in canonical form
 - the index of the row where the non-commuting operator was (that row is now equal to `pauli`; its phase is not updated and for a faithful measurement simulation it needs to be randomized by the user)
 - and the result of the projection if there was no non-commuting operator (`nothing` otherwise)

If `keep_result==false` that result of the projection in case of
anticommutation is not computed, sparing a canonicalization operation. This
canonicalization operation is the only one potentially of cubic complexity.
The rest of the calculations are of quadratic complexity.

If you need to measure a single qubit instead of a multiqubit Pauli operator,
the faster [`projectX!`](@ref), [`projectY!`](@ref), and [`projectZ!`](@ref)
are available.

For less boilerplate and automatic randomization of the phase use [`projectrand!`](@ref).

Here is an example of a projection destroying entanglement:

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

julia> state, anticom_index, result = project!(ghz, P"ZIII");


julia> state
+ Z___
+ Z__Z
+ _Z_Z
+ __ZZ

julia> canonicalize!(state)
+ Z___
+ _Z__
+ __Z_
+ ___Z

julia> anticom_index, result
(1, nothing)
```

And an example of projection consistent with the stabilizer state.

```jldoctest
julia> s = S"ZII
             IXI
             IIY";


julia> canonicalize!(s)
+ _X_
+ __Y
+ Z__

julia> state, anticom_index, result = project!(s, P"-ZII");


julia> state
+ _X_
+ __Y
+ Z__

julia> anticom_index, result
(0, 0x02)
```

While not the best choice, `Stabilizer` can be used for mixed states,
simply by providing an incomplete tableau. In that case it is possible
to attempt to project on an operator that can not be generated by the
provided stabilizer operators. In that case we have `anticom_index==rank`
and `result===nothing`, where `rank` is the the new rank of the tableau,
one more than the number of rows in the initial tableau. However, if
`keep_result` was set to `false`, then `anticom_index` would stay at zero.

```jldoctest
julia> s = S"XZI
             IZI";


julia> project!(s, P"IIX")[1]
+ X__
+ _Z_
```

If we had used [`MixedStabilizer`](@ref) we would have added the projector
to the list of stabilizers.

```jldoctest
julia> s = one(MixedStabilizer, 2, 3)
+ Z__
+ _Z_

julia> project!(s, P"IIX")[1]
+ Z__
+ _Z_
+ __X
```

However, [`MixedDestabilizer`](@ref) would
be an even better choice as it has \$\\mathcal{O}(n^2)\$ complexity
instead of the \$\\mathcal{O}(n^3)\$ complexity of `*Stabilizer`.

```jldoctest
julia> s = one(MixedDestabilizer, 2, 3)
𝒟ℯ𝓈𝓉𝒶𝒷
+ X__
+ _X_
𝒳ₗ━━━
+ __X
𝒮𝓉𝒶𝒷━
+ Z__
+ _Z_
𝒵ₗ━━━
+ __Z

julia> project!(s, P"IIX")[1]
𝒟ℯ𝓈𝓉𝒶𝒷
+ X__
+ _X_
+ __Z
𝒮𝓉𝒶𝒷━
+ Z__
+ _Z_
+ __X
```

See the "Datastructure Choice" section in the documentation for more details.

See also: [`projectX!`](@ref), [`projectY!`](@ref), [`projectZ!`](@ref), [`projectrand!`](@ref)
"""
function project!(state,pauli::PauliOperator;keep_result::Bool=true,phases::Bool=true)
    @valbooldispatch _project!(state,pauli;keep_result=Val(keep_result),phases=Val(phases)) keep_result phases
end

# TODO maybe just add keep_result to it, for consistency
"""
$TYPEDSIGNATURES

When using [`project!`](@ref) on [`MixedStabilizer`](@ref) it automates some of the extra
steps we encounter when implicitly using the `Stabilizer` datastructure to
represent mixed states. Namely, it helps when the projector is not among the
list of stabilizers:

```jldoctest
julia> s = S"XZI
             IZI";


julia> ms = MixedStabilizer(s)
+ X__
+ _Z_

julia> project!(ms, P"IIY")[1]
+ X__
+ _Z_
+ __Y
```

Similarly to [`project!`](@ref) on [`Stabilizer`](@ref), this function has cubic
complexity when the Pauli operator commutes with all rows of the tableau.
Most of the time it is better to simply use [`MixedDestabilizer`](@ref) representation.

Unlike other `project!` methods, this one does not allow for
`keep_result=false`, as the correct rank or anticommutation index can not be
calculated without the expensive (cubic) canonicalization operation required
by `keep_result=true`.

See the "Datastructure Choice" section in the documentation for more details.

See also: [`projectX!`](@ref), [`projectY!`](@ref), [`projectZ!`](@ref).
"""
function project!(state::MixedStabilizer,pauli::PauliOperator;phases::Bool=true)
    @valbooldispatch _project!(state,pauli;phases=Val(phases)) phases
end

function _project!(stabilizer::Stabilizer,pauli::PauliOperator;keep_result::Val{Bkr}=Val(true),phases::Val{Bp}=Val(true)) where {Bkr,Bp}
    anticommutes = 0
    n = size(stabilizer,1)
    for i in 1:n  # The explicit loop is faster than anticommutes = findfirst(row->comm(pauli,stabilizer,row)!=0x0, 1:n); both do not allocate.
        if comm(pauli,stabilizer,i)!=0x0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        if Bkr
            _,_,r = _canonicalize!(stabilizer; phases, ranks=Val(true)) # O(n^3)
            gen = _generate!(copy(pauli), stabilizer; phases) # O(n^2)
            if isnothing(gen)
                result = nothing
                anticommutes = r+1
            else
                result = gen[1].phase[]
            end
        else
            result = nothing
        end
    else
        for i in anticommutes+1:n
            if comm(pauli,stabilizer,i)!=0
                mul_left!(stabilizer, i, anticommutes; phases)
            end
        end
        stabilizer[anticommutes] = pauli
        result = nothing
    end
    stabilizer, anticommutes, result
end

function _project!(ms::MixedStabilizer,pauli::PauliOperator;phases::Val{B}=Val(true)) where B
    _, anticom_index, res = _project!(stabilizerview(ms), pauli; keep_result=Val(true), phases=phases)
    if anticom_index==ms.rank+1 && isnothing(res)
        ms.tab[ms.rank+1] = pauli
        ms.rank += 1
    end
    ms, anticom_index, res
end

@inline function anticomm_update_rows(tab,pauli,r,n,anticommutes,phases::Val{B}=Val(true)) where {B} # TODO Ensure there are no redundant `comm` checks that can be skipped
    for i in r+1:n
        if comm(pauli,tab,i)!=0
            mul_left!(tab, i, n+anticommutes; phases=phases)
        end
    end
    for i in n+anticommutes+1:2n
        if comm(pauli,tab,i)!=0
            mul_left!(tab, i, n+anticommutes; phases=phases)
        end
    end
    for i in 1:r
        if i!=anticommutes && comm(pauli,tab,i)!=0
            mul_left!(tab, i, n+anticommutes; phases=Val(false))
        end
    end
end

function _project!(d::Destabilizer,pauli::PauliOperator;keep_result::Val{Bkr}=Val(true),phases::Val{Bp}=Val(true)) where {Bkr, Bp} # repetition between Destabilizer and MixedDestabilizer, but the redundancy makes the two codes slightly simpler and easier to infer
    anticommutes = 0
    tab = d.tab
    stabilizer = stabilizerview(d)
    destabilizer = destabilizerview(d)
    r = trusted_rank(d)
    n = length(d) # not `nqubits(d)` in case we have an incomplete tableau    # Check whether we anticommute with any of the stabilizer rows
    for i in 1:r # The explicit loop is faster than anticommutes = findfirst(row->comm(pauli,stabilizer,row)!=0x0, 1:r); both do not allocate.
        if comm(pauli,stabilizer,i)!=0x0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        if n != nqubits(stabilizer)
            throw(BadDataStructure("`Destabilizer` can not efficiently (faster than n^3) detect whether you are projecting on a stabilized or a logical operator. Switch to one of the `Mixed*` data structures.",
                                   :project!,
                                   :Destabilizer))
        end
        if Bkr
            new_pauli = zero(pauli)::PauliOperator # typeassert for JET
            new_pauli.phase[] = pauli.phase[]
            for i in 1:r
                comm(pauli,destabilizer,i)!=0 && mul_left!(new_pauli, stabilizer, i, phases=phases)
            end
            result = new_pauli.phase[]
        else
            result = nothing
        end
    else
        anticomm_update_rows(tab,pauli,r,n,anticommutes,phases)
        destabilizer[anticommutes] = stabilizer[anticommutes]
        stabilizer[anticommutes] = pauli
        result = nothing
    end
    d, anticommutes, result
end

function _project!(d::MixedDestabilizer,pauli::PauliOperator;keep_result::Val{Bkr}=Val(true),phases::Val{Bp}=Val(true)) where {Bkr, Bp} # repetition between Destabilizer and MixedDestabilizer, but the redundancy makes the two codes slightly simpler and easier to infer
    anticommutes = 0
    tab = d.tab
    stabilizer = stabilizerview(d)
    destabilizer = destabilizerview(d)
    r = trusted_rank(d)
    n = length(d) # not `nqubits(d)` in case we have an incomplete tableau    # Check whether we anticommute with any of the stabilizer rows
    for i in 1:r # The explicit loop is faster than anticommutes = findfirst(row->comm(pauli,stabilizer,row)!=0x0, 1:r); both do not allocate.
        if comm(pauli,stabilizer,i)!=0x0
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        anticomlog = 0
        # Check whether we anticommute with any of the logical X rows
        for i in r+1:n # The explicit loop is faster than findfirst.
            if comm(pauli,tab,i)!=0x0
                anticomlog = i
                break
            end
        end
        if anticomlog==0
            # Check whether we anticommute with any of the logical Z rows
            for i in n+r+1:2*n # The explicit loop is faster than findfirst.
                if comm(pauli,tab,i)!=0x0
                    anticomlog = i
                    break
                end
            end
        end
        if anticomlog!=0
            if anticomlog<=n
                rowswap!(tab, r+1+n, anticomlog)
                n!=r+1 && anticomlog!=r+1 && rowswap!(tab, r+1, anticomlog+n)
            else
                rowswap!(tab, r+1, anticomlog-n)
                rowswap!(tab, r+1+n, anticomlog)
            end
            anticomm_update_rows(tab,pauli,r+1,n,r+1,phases)
            d.rank+=1
            anticommutes = d.rank
            tab[r+1] = tab[n+r+1]
            tab[n+r+1] = pauli
            result = nothing
        else
            if Bkr
                new_pauli = zero(pauli)::PauliOperator # typeassert for JET
                new_pauli.phase[] = pauli.phase[]
                for i in 1:r
                    comm(pauli,destabilizer,i)!=0 && mul_left!(new_pauli, stabilizer, i, phases=phases)
                end
                result = new_pauli.phase[]
            else
                result = nothing
            end
        end
    else
        anticomm_update_rows(tab,pauli,r,n,anticommutes,phases)
        destabilizer[anticommutes] = stabilizer[anticommutes]
        stabilizer[anticommutes] = pauli
        result = nothing
    end
    d, anticommutes, result
end

"""
Measure a given qubit in the X basis.
A faster special-case version of [`project!`](@ref).

See also: [`project!`](@ref), [`projectXrand!`](@ref), [`projectY!`](@ref), [`projectZ!`](@ref).
"""
function projectX!(d::MixedDestabilizer,qubit::Int;keep_result::Bool=true,phases::Bool=true)
    @valbooldispatch project_cond!(d,qubit,Val(isZ),Val((true,false));keep_result,phases=Val(phases)) phases
end

"""
Measure a given qubit in the Z basis.
A faster special-case version of [`project!`](@ref).

See also: [`project!`](@ref), [`projectZrand!`](@ref), [`projectY!`](@ref), [`projectX!`](@ref).
"""
function projectZ!(d::MixedDestabilizer,qubit::Int;keep_result::Bool=true,phases::Bool=true)
    @valbooldispatch project_cond!(d,qubit,Val(isX),Val((false,true));keep_result,phases=Val(phases))
end

"""
Measure a given qubit in the Y basis.
A faster special-case version of [`project!`](@ref).

See also: [`project!`](@ref), [`projectYrand!`](@ref), [`projectX!`](@ref), [`projectZ!`](@ref).
"""
function projectY!(d::MixedDestabilizer,qubit::Int;keep_result::Bool=true,phases::Bool=true)
    @valbooldispatch project_cond!(d,qubit,Val(isXorZ),Val((true,true));keep_result,phases=Val(phases))
end

@inline isX(tab,row,col) = tab[row,col][1]
@inline isZ(tab,row,col) = tab[row,col][2]
@inline isY(tab,row,col) = (&)(tab[row,col]...)
@inline isXorZ(tab,row,col) = ⊻(tab[row,col]...)

@inline function anticomm_update_rows_cond(tab,qubit,r,n,anticommutes,phases::Val{B},cond::Val{IS}) where {B,IS} # TODO Ensure there are no redundant `comm` checks that can be skipped
    for i in r+1:n
        if IS(tab,i,qubit)
            mul_left!(tab, i, n+anticommutes; phases=phases)
        end
    end
    for i in n+anticommutes+1:2n
        if IS(tab,i,qubit)
            mul_left!(tab, i, n+anticommutes; phases=phases)
        end
    end
    for i in 1:r
        if i!=anticommutes && IS(tab,i,qubit)
            mul_left!(tab, i, n+anticommutes; phases=Val(false))
        end
    end
end

"""Internal method used to implement [`projectX!`](@ref), [`projectZ!`](@ref), and [`projectY!`](@ref)."""
function project_cond!(d::MixedDestabilizer,qubit::Int,cond::Val{IS},reset::Val{RESET};keep_result::Bool=true,phases::Val{PHASES}=Val(true)) where {IS,RESET,PHASES}
    anticommutes = 0
    tab = d.tab
    stabilizer = stabilizerview(d)
    destabilizer = destabilizerview(d)
    r = d.rank
    n = nqubits(d)
    # Check whether we anticommute with any of the stabilizer rows
    for i in 1:r # The explicit loop is faster than anticommutes = findfirst(row->comm(pauli,stabilizer,row)!=0x0, 1:r); both do not allocate.
        if IS(stabilizer,i,qubit)
            anticommutes = i
            break
        end
    end
    if anticommutes == 0
        anticomlog = 0
        # Check whether we anticommute with any of the logical X rows
        for i in r+1:n # The explicit loop is faster than findfirst.
            if IS(tab,i,qubit)
                anticomlog = i
                break
            end
        end
        if anticomlog==0
            # Check whether we anticommute with any of the logical Z rows
            for i in n+r+1:2*n # The explicit loop is faster than findfirst.
                if IS(tab,i,qubit)
                    anticomlog = i
                    break
                end
            end
        end
        if anticomlog!=0
            if anticomlog<=n
                rowswap!(tab, r+1+n, anticomlog)
                n!=r+1 && anticomlog!=r+1 && rowswap!(tab, r+1, anticomlog+n)
            else
                rowswap!(tab, r+1, anticomlog-n)
                rowswap!(tab, r+1+n, anticomlog)
            end
            anticomm_update_rows_cond(tab,qubit,r+1,n,r+1,phases,cond)
            d.rank += 1
            anticommutes = d.rank
            tab[r+1] = tab[n+r+1]
            zero!(tab,n+r+1) # set to projector
            tab[n+r+1,qubit] = RESET
            result = nothing
        else
            if keep_result
                new_pauli = zero(PauliOperator,n)
                for i in 1:r # comm check bellow
                    IS(destabilizer,i,qubit) && mul_left!(new_pauli, stabilizer, i, phases=phases)
                end
                result = new_pauli.phase[]
            else
                result = nothing
            end
        end
    else
        anticomm_update_rows_cond(tab,qubit,r,n,anticommutes,phases,cond)
        destabilizer[anticommutes] = stabilizer[anticommutes]
        zero!(stabilizer, anticommutes) # set to projector
        stabilizer[anticommutes,qubit] = RESET
        result = nothing
    end
    d, anticommutes, result
end

"""
$TYPEDSIGNATURES

Project `qubit` of `state` along the X axis and randomize the phase if necessary.

Lower boilerplate version of [`project!`](@ref).

See also: [`project!`](@ref), [`projectX!`](@ref), [`projectZrand!`](@ref), [`projectYrand!`](@ref)
"""
function projectXrand!(state, qubit)
    _, anticom, res = projectX!(state, qubit)
    isnothing(res) && (res = tab(stabilizerview(state)).phases[anticom] = rand((0x0, 0x2)))
    return state, res
end

"""
$TYPEDSIGNATURES

Project `qubit` of `state` along the Z axis and randomize the phase if necessary.

Lower boilerplate version of [`project!`](@ref).

See also: [`project!`](@ref), [`projectZ!`](@ref), [`projectXrand!`](@ref), [`projectYrand!`](@ref)
"""
function projectZrand!(state, qubit)
    _, anticom, res = projectZ!(state, qubit)
    isnothing(res) && (res = tab(stabilizerview(state)).phases[anticom] = rand((0x0, 0x2)))
    return state, res
end

"""
$TYPEDSIGNATURES

Project `qubit` of `state` along the Y axis and randomize the phase if necessary.

Lower boilerplate version of [`project!`](@ref).

See also: [`project!`](@ref), [`projectY!`](@ref), [`projectXrand!`](@ref), [`projectZrand!`](@ref)
"""
function projectYrand!(state, qubit)
    _, anticom, res = projectY!(state, qubit)
    isnothing(res) && (res = tab(stabilizerview(state)).phases[anticom] = rand((0x0, 0x2)))
    return state, res
end

"""
$TYPEDSIGNATURES

Measure `pauli` operator on `state` and randomize the phase if necessary.

Lower boilerplate version of [`project!`](@ref).

See also: [`project!`](@ref), [`projectXrand!`](@ref), [`projectZrand!`](@ref), [`projectYrand!`](@ref)
"""
function projectrand!(state, pauli)
    _, anticom, res = project!(state, pauli)
    isnothing(res) && (res = tab(stabilizerview(state)).phases[anticom] = rand((0x0, 0x2)))
    return state, res
end

"""
$TYPEDSIGNATURES

Trace out a qubit.

See also: [`delete_columns`](@ref)

""" # TODO all of these should raise an error if length(qubits)>rank
function traceout!(s::Stabilizer, qubits; phases=true, rank=false)
    _,i = canonicalize_rref!(s,qubits;phases=phases)
    for j in i+1:size(s,1)
        zero!(s,j)
    end
    if rank return (s, i) else return s end
end

"""
$TYPEDSIGNATURES
"""
function traceout!(s::Union{MixedStabilizer, MixedDestabilizer}, qubits; phases=true, rank=false)
    _,i = canonicalize_rref!(s,qubits;phases=phases)
    s.rank = i
    if rank return (s, i) else return s end
end

function _expand_pauli(pauli,qubits,n) # TODO rename and make public
    expanded = zero(PauliOperator,n)
    for (ii, i) in enumerate(qubits)
        expanded[i] = pauli[ii]
    end
    expanded.phase[] = pauli.phase[]
    expanded
end

"""
$TYPEDSIGNATURES

Reset a given set of qubits to be in the state `newstate`.
These qubits are traced out first, which could lead to "nonlocal" changes in
the tableau.
"""
function reset_qubits!(s::Stabilizer, newstate, qubits; phases=true)
    nqubits(newstate)==length(qubits) || throw(DimensionMismatch("`qubits` and `newstate` have to be of consistent size"))
    length(qubits) <= nqubits(s) || throw(DimensionMismatch("the stabilizer is not big enough to contain the new state"))
    n = nqubits(s)
    s, x, z = canonicalize!(s,ranks=true) # TODO this is unnecessary, but it provides for nicely formatted tableaux; consider removing it for speed reasons
    _, rref_i = canonicalize_rref!((@view s[1:z]),qubits,phases=phases)
    for row in 1:length(newstate)
        s[row+rref_i] = _expand_pauli(newstate[row], qubits, n) # TODO do something that does not allocate temporary arrays
    end
    for row in rref_i+length(newstate)+1:z
        zero!(s,row)
    end
    s
end

"""
$TYPEDSIGNATURES
"""
function reset_qubits!(s::MixedStabilizer, newstate, qubits; phases=true) # TODO create the necessary interfaces so that Stabilizer and MixedStabilizer share this code
    nqubits(newstate)==length(qubits) || throw(DimensionMismatch("`qubits` and `newstate` have to be of consistent size"))
    length(qubits) <= nqubits(s) || throw(DimensionMismatch("the stabilizer is not big enough to contain the new state"))
    n = nqubits(s)
    sv = stabilizerview(s)
    sv, rref_i = canonicalize_rref!(sv,qubits,phases=phases)
    for row in 1:length(newstate)
        s.tab[row+rref_i] = _expand_pauli(newstate[row], qubits, n) # TODO do something that does not allocate temporary arrays
    end
    s.rank = rref_i+length(newstate)
    s
end

"""
$TYPEDSIGNATURES
"""
function reset_qubits!(s::MixedDestabilizer, newstate::AbstractStabilizer, qubits; phases=true) # TODO this is really inefficient
    _phases = Val(phases)
    nqubits(newstate)==length(qubits) || throw(DimensionMismatch("`qubits` and `newstate` have to be of consistent size"))
    length(qubits) <= nqubits(s) || throw(DimensionMismatch("the stabilizer is not big enough to contain the new state"))
    newstatestab = stabilizerview(newstate)
    traceout!(s,qubits)
    for pauli in newstatestab
        expanded = _expand_pauli(pauli, qubits, nqubits(s)) # TODO, use a sparse project that does not require this expand
        _, anticomm, res = project!(s,expanded, phases=phases) # TODO make an `apply_measurement_phase!(project!(...), phase)`
        sv =  stabilizerview(s)
        if anticomm!=0 # Does not commute with the stabilizer or logical ops
            tab(sv).phases[anticomm] = pauli.phase[]
        else # Commutes with everyone
            if res!=0 && phases # TODO many of the checks below were already done by project!; find a way to not repeat them
                destab = destabilizerview(s)
                r = LinearAlgebra.rank(s)
                loc = findfirst(i->comm(pauli,destab,i)!=0, 1:r)::Int # `nothing` should not be a possible answer
                for i in loc+1:r
                    if comm(pauli, destab, i)!=0
                        mul_left!(s, i, loc; phases=_phases)
                    end
                end
                sv[loc] = pauli
            end
        end
    end
    s
end

"""
    expect(p::PauliOperator, st::AbstractStabilizer)

Compute the expectation value of a Pauli operator `p` on a stabilizer state `st`.
This function will allocate a temporary copy of the stabilizer state `st`.
"""
function expect(p::PauliOperator, s::AbstractStabilizer)
    nqubits(p) == nqubits(s) || error("The number of qubits does not match")
    _, _, result = project!(copy(s), p)
    result === nothing && return 0
    result === 0x00 && return 1
    result === 0x01 && return im
    result === 0x02 && return -1
    result === 0x03 && return -im
end

"""Put source tableau in target tableau at given row and column. Assumes target location is zeroed out.""" # TODO implement a getindex setindex interface to this
@inline function puttableau!(target::Tableau{V1,M1}, source::Tableau{V2,M2}, row::Int, col::Int; phases::Val{B}=Val(true)) where {B,V1,V2,T<:Unsigned,M1<:AbstractMatrix{T},M2<:AbstractMatrix{T}}
    xzs = target.xzs
    ph = target.phases
    sxzs = source.xzs
    sph = source.phases
    r,n = size(source)
    bₗ = _div(T,col) + 1
    bᵣ = bₗ + 1
    eₗ = bₗ + _div(T,n-1)
    eᵣ = _div(T,col+n-1) + 1
    shiftₗ = _mod(T,col)
    shiftᵣ = 8*sizeof(T)-shiftₗ
    for i in 1:r
    @inbounds @simd for j in 0:eₗ-bₗ # TODO more simdification
        xzs[bₗ+j,row+i] |= sxzs[j+1,i] >>> -shiftₗ
        xzs[end÷2+bₗ+j,row+i] |= sxzs[end÷2+j+1,i] >>> -shiftₗ
    end
    @inbounds @simd for j in 0:eᵣ-bᵣ
        xzs[bᵣ+j,row+i] |= sxzs[j+1,i] >>> shiftᵣ
        xzs[end÷2+bᵣ+j,row+i] |= sxzs[end÷2+j+1,i] >>> shiftᵣ
    end
    end
    B && (ph[row+1:row+r] .= sph)
    target, row+r, col+n
end

function puttableau!(target::Stabilizer{T1}, source::Stabilizer{T2}, row::Int, col::Int; phases::Val{B}=Val(true)) where {B,T1,T2}
    puttableau!(tab(target), tab(source), row, col; phases)
end
function puttableau!(target::Tableau, source::Stabilizer{T2}, row::Int, col::Int; phases::Val{B}=Val(true)) where {B,T2}
    puttableau!(target, tab(source), row, col; phases)
end
function puttableau!(target::Stabilizer{T1}, source::Tableau, row::Int, col::Int; phases::Val{B}=Val(true)) where {B,T1}
    puttableau!(tab(target), source, row, col; phases)
end

"""Unexported low-level function that removes a column (by shifting all columns to the right of the target by one step to the left)

Because Tableau is not mutable we return a new Tableau with the same (modified) xzs array."""
function remove_column!(s::Tableau{V,M}, col::Int) where {V,T<:Unsigned,M<:AbstractMatrix{T}}
    rows,cols=size(s)
    xzs = s.xzs
    big = _div(T,col-1) + 1
    nbig = size(s.xzs,1)÷2
    shiftᵣ = 8*sizeof(T)-1
    zero_span = ~zero(T) << _mod(T,col-1)
    zero_first = ~zero(T) >> 1
    for i in 1:rows # TODO more simdification
        xzs[big,i]        = (xzs[big,i]       & ~zero_span) | ( zero_span & (xzs[big,i]       >>> 1))
        xzs[end÷2+big,i]  = (xzs[end÷2+big,i] & ~zero_span) | ( zero_span & (xzs[end÷2+big,i] >>> 1))
        for j in big+1:nbig
            xzs[j-1,i]        = (xzs[j-1,i]       & zero_first) | (xzs[j,i]       << shiftᵣ)
            xzs[end÷2+j-1,i]  = (xzs[end÷2+j-1,i] & zero_first) | (xzs[end÷2+j,i] << shiftᵣ)
            xzs[j,i] = xzs[j,i] >> 1
            xzs[end÷2+j,i] = xzs[end÷2+j,i] >> 1
        end
    end
    Tableau(s.phases, nqubits(s)-1, s.xzs)
end
remove_column!(s::Stabilizer, col::Int) = Stabilizer(remove_column!(tab(s),col))

"""Unexported low-level function that moves row i to row j.

Used on its own, this function will break invariants. Meant to be used in `_remove_rowcol!`.
"""
@inline function _rowmove!(s::Tableau, i, j; phases::Val{B}=Val(true)) where B
    (i == j) && return
    B && begin s.phases[j] = s.phases[i] end
    @inbounds @simd for k in 1:size(s.xzs,1)
        s.xzs[k,j] = s.xzs[k,i]
    end
end
@inline _rowmove!(s::Stabilizer, i, j; phases::Val{B}=Val(true)) where B = _rowmove!(tab(s), i, j; phases)

"""Unexported low-level function that removes a row (by shifting all rows up as necessary)

Because MixedDestabilizer is not mutable we return a new MixedDestabilizer with the same (modified) xzs array.

Used on its own, this function will break invariants. Meant to be used with `projectremove!`.
"""
function _remove_rowcol!(s::MixedDestabilizer, r,c)
    t = tab(s)
    rows, cols = size(t)
    t = remove_column!(t,c) # TODO col and row removal should be done in the same loop, not two separate loops
    for i in r+1:cols+r-1
        _rowmove!(t, i, i-1)
    end
    for i in cols+r+1:rows
        _rowmove!(t, i, i-2)
    end
    oldrank = rank(s)
    newrank = r<=oldrank ? oldrank-1 : oldrank
    MixedDestabilizer(Tableau((@view t.phases[1:end-2]),cols-1,(@view t.xzs[:,1:end-2])), newrank)
end

#=
"""Unexported low-level function that projects a qubit and returns the result while making the tableau smaller by a qubit.

Because MixedDestabilizer is not mutable we return a new MixedDestabilizer with the same (modified) xzs array.
"""
function projectremove!(s::MixedDestabilizer, projfunc::F, qubit) where {F<:Union{typeof(projectX!),typeof(projectY!),typeof(projectZ!)}}
    error("can not be implemented in the style of project!, because one can not change `res` after the row has been removed")
end
=#

"""Unexported low-level function that projects a qubit and returns the result while making the tableau smaller by a qubit.

Because MixedDestabilizer is not mutable we return a new MixedDestabilizer with the same (modified) xzs array.
"""
function projectremoverand!(s::MixedDestabilizer, projfunc::F, qubit) where {F<:Union{typeof(projectX!),typeof(projectY!),typeof(projectZ!)}}
    _, anticom, res = projfunc(s, qubit)
    if anticom!=0
        res = rand((0x0,0x2))
        phases(stabilizerview(s))[anticom] = res
    end
    r = rank(s)
    traceout!(s,[qubit]) # TODO this can be optimized thanks to the information already known from projfunc
    s = _remove_rowcol!(s, r, qubit)
    s, res
end

function traceoutremove!(s::MixedDestabilizer, qubit)
    traceout!(s,[qubit]) # TODO this can be optimized thanks to the information already known from projfunc
    s = _remove_rowcol!(s, nqubits(s), qubit)
end
