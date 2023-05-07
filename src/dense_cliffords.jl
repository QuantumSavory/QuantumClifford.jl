"""
Clifford Operator specified by the mapping of the basis generators.

```jldoctest
julia> tCNOT
X₁ ⟼ + XX
X₂ ⟼ + _X
Z₁ ⟼ + Z_
Z₂ ⟼ + ZZ

julia> phase_gate = C"Y
                      Z"
X₁ ⟼ + Y
Z₁ ⟼ + Z

julia> stab = S"XI
                IZ";


julia> entangled = tCNOT*stab
+ XX
+ ZZ
```

julia> CliffordOperator(T"YY")
ERROR: DimensionMismatch("Input tableau should be square (in which case the destabilizers are calculated) or of size 2n×n (in which case it is used directly).")
[...]
```

[`Destabilizer`](@ref) can also be converted.
```jldoctest
julia> d = Destabilizer(S"Y")
+ Z
━━━
+ Y

julia> CliffordOperator(d)
X ⟼ + Z
Z ⟼ + Y
```
"""
struct CliffordOperator{T<:Tableau} <: AbstractCliffordOperator
    tab::T
    function CliffordOperator(tab::Tableau)
        if size(tab,1)==2*size(tab,2)
            new{typeof(tab)}(tab)
        #elseif size(stab,1)==size(stab,2) # TODO be able to work with squara tableaux (by reversing all row operations)
        #    destab = tab(Destabilizer(stab))
        #    new{typeof(destab.phases),typeof(destab.xzs)}(destab) # TODO be smarter about type signatures here... there should be a better way
        else
            throw(DimensionMismatch("Input tableau should be of size 2n×n (top half is the X mappings and the bottom half are the Z mappings)."))
        end
    end
end

macro C_str(a)
    tab = _T_str(a)
    quote CliffordOperator($tab) end
end

CliffordOperator(op::CliffordOperator) = op
CliffordOperator(paulis::AbstractVector{<:PauliOperator}) = CliffordOperator(Tableau(paulis))
CliffordOperator(destab::Destabilizer) = CliffordOperator(tab(destab))

function CliffordOperator(stab::Stabilizer)
    md = MixedDestabilizer(copy(stab); trackoperations=true)
    cltab = tab(md)
    # cltab = vcat(cltab[1+length(cltab)÷2:length(cltab)], cltab[1:length(cltab)÷2])
    return CliffordOperator(cltab)
end

Base.:(==)(l::CliffordOperator, r::CliffordOperator) = l.tab == r.tab
Base.hash(c::T, h::UInt) where {T<:CliffordOperator} = hash(T, hash(tab(c), h))

Base.getindex(c::CliffordOperator, args...) = getindex(tab(c), args...)
Base.setindex!(c::CliffordOperator, args...) = setindex!(tab(c), args...)

tab(c::CliffordOperator) = c.tab

Base.size(c::CliffordOperator,args...) = size(tab(c),args...)

function row_limit(str, limit=50)
    n = length(str)
    if (n <= limit || limit==-1)
        return str
    end
    padding = Int64(floor(limit/2))
    return SubString(str, 1, padding) * " ... " * SubString(str, n - padding, n - 1)
end

digits_subchars = collect("₀₁₂₃₄₅₆₇₈₉")
digits_substr(n,nwidth) = join(([digits_subchars[d+1] for d in reverse(digits(n, pad=nwidth))]))

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

function Base.copy(c::CliffordOperator)
    CliffordOperator(copy(c.tab))
end

@inline nqubits(c::CliffordOperator) = nqubits(c.tab)

Base.zero(c::CliffordOperator) = CliffordOperator(zero(c.tab))
Base.zero(::Type{<:CliffordOperator}, n) = CliffordOperator(zero(Tableau, 2n, n))

function Base.:(*)(l::AbstractCliffordOperator, r::CliffordOperator)
    tab = copy(r.tab)
    apply!(Stabilizer(tab),l) # TODO maybe not the most elegant way to perform apply!(::Tableau, gate)
    CliffordOperator(tab)
end

function apply!(r::CliffordOperator, l::AbstractCliffordOperator; phases=false)
    @valbooldispatch _apply!(Stabilizer(tab(r)),l;phases=Val(phases)) phases # TODO maybe not the most elegant way to perform apply!(::Tableau, gate)
    r
end

# TODO create Base.permute! and getindex(..., permutation_array)
function permute(c::CliffordOperator,p) # TODO this is a slow stupid implementation
    CliffordOperator(Tableau([c.tab[i][p] for i in 1:2*nqubits(c)][vcat(p,p.+nqubits(c))]))
end

"""Nonvectorized version of `apply!` used for unit tests."""
function _apply_nonthread!(stab::AbstractStabilizer, c::CliffordOperator; phases::Bool=true)
    nqubits(stab)==nqubits(c) || throw(DimensionMismatch("The tableau and the Clifford operator need to act on the same number of qubits. Consider specifying an array of indices as a third argument to the `apply!` function to avoid this error."))
    s_tab = tab(stab)
    c_tab = tab(c)
    new_stabrow = zero(s_tab[1])
    for row_stab in eachindex(s_tab)
        zero!(new_stabrow)
        apply_row_kernel!(new_stabrow, row_stab, s_tab, c_tab, phases=Val(phases))
    end
    stab
end

# TODO no need to track phases outside of stabview
function _apply!(stab::AbstractStabilizer, c::CliffordOperator; phases::Val{B}=Val(true)) where B
    nqubits(stab)==nqubits(c) || throw(DimensionMismatch("The tableau and the Clifford operator need to act on the same number of qubits. Consider specifying an array of indices as a third argument to the `apply!` function to avoid this error."))
    s_tab = tab(stab)
    c_tab = tab(c)
    @batch minbatch=25 threadlocal=zero(c_tab[1]) for row_stab in eachindex(s_tab)
        zero!(threadlocal) # a new stabrow for temporary storage
        apply_row_kernel!(threadlocal, row_stab, s_tab, c_tab, phases=phases)
    end
    stab
end

# TODO Added a lot of type assertions to help Julia infer types, but they are much too strict for cases where bitpacking varies (check tests)
#@inline function apply_row_kernel!(new_stabrow::PauliOperator{Array{UInt8,0},Vector{Tme}}, row::Int, s_tab::Stabilizer{Tv,Tm}, c_tab::Stabilizer{Tv,Tm}; phases=true) where {Tme,Tv<:AbstractVector{UInt8},Tm<:AbstractMatrix{Tme}}
@inline function apply_row_kernel!(new_stabrow, row, s_tab, c_tab; phases::Val{B}=Val(true)) where B
    B && (new_stabrow.phase[] = s_tab.phases[row])
    n = nqubits(c_tab)
    for qubit in 1:n
        x,z = s_tab[row,qubit]
        if B&&x&&z
            new_stabrow.phase[] -= 0x1
        end
        if x
            mul_left!(new_stabrow, c_tab, qubit, phases=phases)
        end
        if z
            mul_left!(new_stabrow, c_tab, qubit+n, phases=phases)
        end
    end
    s_tab[row] = new_stabrow
    new_stabrow
end

"""Nonvectorized version of `apply!` used for unit tests."""
function _apply_nonthread!(stab::AbstractStabilizer, c::CliffordOperator, indices_of_application::AbstractArray{Int,1}; phases::Bool=true)
    s_tab = tab(stab)
    c_tab = tab(c)
    new_stabrow = zero(PauliOperator,nqubits(c))
    for row in eachindex(s_tab)
        zero!(new_stabrow)
        apply_row_kernel!(new_stabrow, row, s_tab, c_tab, indices_of_application; phases=Val(phases))
    end
    stab
end

#TODO a lot of code repetition with apply!(stab::AbstractStabilizer, c::CliffordOperator; phases::Bool=true) and apply_row_kernel!
function _apply!(stab::AbstractStabilizer, c::CliffordOperator, indices_of_application::AbstractArray{Int,1}; phases::Val{B}=Val(true)) where B
    #max(indices_of_application)<=nqubits(s) || throw(DimensionMismatch("")) # Too expensive to check every time
    s_tab = tab(stab)
    c_tab = tab(c)
    @batch minbatch=25 threadlocal=zero(c_tab[1]) for row_stab in eachindex(s_tab)
        zero!(threadlocal) # a new stabrow for temporary storage
        apply_row_kernel!(threadlocal, row_stab, s_tab, c_tab, indices_of_application, phases=phases)
    end
    stab
end

@inline function apply_row_kernel!(new_stabrow, row, s_tab, c_tab, indices_of_application; phases::Val{B}=Val(true)) where B
    B && (new_stabrow.phase[] = s_tab.phases[row])
    n = nqubits(c_tab)
    for (qubit_i, qubit) in enumerate(indices_of_application)
        x,z = s_tab[row,qubit]
        if B&&x&&z
            new_stabrow.phase[] -= 0x1
        end
        if x
            mul_left!(new_stabrow, c_tab, qubit_i, phases=phases)
        end
        if z
            mul_left!(new_stabrow, c_tab, qubit_i+n, phases=phases)
        end
    end
    for (qubit_i, qubit) in enumerate(indices_of_application)
        s_tab[row,qubit] = new_stabrow[qubit_i]
    end
    B && (s_tab.phases[row] = new_stabrow.phase[])
    new_stabrow
end

const tCNOT = C"XX
                IX
                ZI
                ZZ"

const tCPHASE = C"XZ
                  ZX
                  ZI
                  IZ"

const tSWAP = C"IX
                XI
                IZ
                ZI"

const tHadamard = C"Z
                    X"

const tPhase = C"Y
                 Z"

const tId1 = C"X
               Z"
