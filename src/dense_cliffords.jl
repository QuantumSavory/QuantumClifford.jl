"""
Clifford Operator specified by the mapping of the basis generators.

```jldoctest
julia> CNOT
X_ ⟼ + XX
_X ⟼ + _X
Z_ ⟼ + Z_
_Z ⟼ + ZZ

julia> phase_gate = C"Y
                      Z"
X ⟼ + Y
Z ⟼ + Z

julia> stab = S"XI
                IZ";

julia> entangled = CNOT*stab
+ XX
+ ZZ
```

You can convert a stabilizer into a Clifford Operator (if necessary, the destabilizers are calculated on the fly):

```jldoctest
julia> CliffordOperator(S"Y")
X ⟼ + Z
Z ⟼ + Y


julia> CliffordOperator(S"YY")
ERROR: DimensionMismatch("Input tableau should be square (in which case the destabilizers are calculated) or of size 2n×n (in which case it is used directly).")
[...]
```

[`Destabilizer`](@ref) can also be converted (actually, internally, square stabilizer tableaux are first converted to destabilizer tableaux).
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
struct CliffordOperator{Tzv<:AbstractVector{UInt8},Tm<:AbstractMatrix{<:Unsigned}} <: AbstractCliffordOperator
    tab::Stabilizer{Tzv,Tm}
    function CliffordOperator(stab::Stabilizer{Tzv,Tm}) where {Tzv,Tm}
        if size(stab,1)==2*size(stab,2)
            new{Tzv,Tm}(stab)
        #elseif size(stab,1)==size(stab,2) # TODO be able to work with squara tableaux (by reversing all row operations)
        #    destab = tab(Destabilizer(stab))
        #    new{typeof(destab.phases),typeof(destab.xzs)}(destab) # TODO be smarter about type signatures here... there should be a better way
        else
            throw(DimensionMismatch("Input tableau should be of size 2n×n (top half is the X mappings and the bottom half are the Z mappings)."))
        end
    end
end

macro C_str(a)
    tab = _S_str(a)
    CliffordOperator(tab)
end

CliffordOperator(op::CliffordOperator) = op
CliffordOperator(paulis::AbstractVector{<:PauliOperator}) = CliffordOperator(Stabilizer(paulis))
CliffordOperator(destab::Destabilizer) = CliffordOperator(tab(destab))

Base.:(==)(l::CliffordOperator, r::CliffordOperator) = l.tab == r.tab
Base.hash(c::T, h::UInt) where {T<:CliffordOperator} = hash(T, hash(tab(c), h))

Base.getindex(c::CliffordOperator, args...) = getindex(tab(c), args...)
Base.setindex!(c::CliffordOperator, args...) = setindex!(tab(c), args...)

tab(c::CliffordOperator) = c.tab

Base.size(c::CliffordOperator,args...) = size(tab(c),args...)

function Base.show(io::IO, c::CliffordOperator)
    n = nqubits(c)
    for i in 1:n
        print(io, repeat("_",i-1),"X",repeat("_",n-i), " ⟼ ")
        print(io, c.tab[i])
        println(io)
    end
    for i in 1:n
        print(io, repeat("_",i-1),"Z",repeat("_",n-i), " ⟼ ")
        print(io, c.tab[i+n])
        println(io)
    end
end

function Base.copy(c::CliffordOperator)
    CliffordOperator(copy(c.tab))
end

@inline nqubits(c::CliffordOperator) = nqubits(c.tab)

Base.zero(c::CliffordOperator) = CliffordOperator(zero(c.tab))
Base.zero(::Type{<:CliffordOperator}, n) = CliffordOperator(zero(Stabilizer, 2n, n))

function Base.:(*)(l::AbstractCliffordOperator, r::CliffordOperator)
    tab = copy(r.tab)
    apply!(tab,l)
    CliffordOperator(tab)
end

function apply!(r::CliffordOperator, l::AbstractCliffordOperator; phases=false)
    _apply!(tab(r),l;phases=Val(phases))
    r
end

# TODO create Base.permute! and getindex(..., permutation_array)
function permute(c::CliffordOperator,p::AbstractArray{T,1} where T) # TODO this is a slow stupid implementation
    CliffordOperator(Stabilizer([c.tab[i][p] for i in 1:2*nqubits(c)][vcat(p,p.+nqubits(c))]))
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
