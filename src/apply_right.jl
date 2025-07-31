"""
the `apply_right!` function is used to right multiply any quantum operation to unitary 
Clifford operation or Pauli product

```jldoctest
julia> apply_right!(C"X Z", sHadamard(1))
X₁ ⟼ + Z
Z₁ ⟼ + X
julia> apply_right!(C"Y Z", C"Z Y")
X₁ ⟼ + Z
Z₁ ⟼ - X
julia> apply_right!(C"Y Z", P"X")
X₁ ⟼ + Y
Z₁ ⟼ - Z
```

Example: Build a bell state decoder
```jldoctest
julia> cliff = one(CliffordOperator, 2)
X₁ ⟼ + X_
X₂ ⟼ + _X
Z₁ ⟼ + Z_
Z₂ ⟼ + _Z
julia> apply_right!(cliff, sHadamard(1))
X₁ ⟼ + Z_
X₂ ⟼ + _X
Z₁ ⟼ + X_
Z₂ ⟼ + _Z
julia> apply_right!(cliff, sCNOT(1, 2))
X₁ ⟼ + ZX
X₂ ⟼ + _X
Z₁ ⟼ + X_
Z₂ ⟼ + XZ
julia> apply!(bell(), cliff)
+ Z_
+ _Z
```

See also: [`apply!`](@ref), [`apply_inv!`](@ref)
"""
function apply_right! end


##############################
# Dense operators
##############################

function apply_right!(l::CliffordOperator, r::PauliOperator)
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The Clifford and Pauli operators need to act on the same number of qubits."))
    for i in 1:nqubits(l)
        x, z = r[i]
        if x
            phases(l)[nqubits(l)+i] ⊻= 0x02
        end
        if z
            phases(l)[i] ⊻= 0x02
        end
    end
    return l
end

function apply_right!(l::CliffordOperator, r::CliffordOperator; phases=true)
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The two Clifford operators need to act on the same number of qubits."))
    l_tab = tab(l)
    l_tab_copy = copy(l_tab)
    r_tab = tab(r)
    threadlocal = l.buffer
    @inbounds for row_r in eachindex(r_tab)
        zero!(threadlocal)
        apply_right_row_kernel!(threadlocal, l_tab, row_r, l_tab_copy, r_tab; phases)
    end
    l
end

"""helper for computing the right multiplication of a row of a Clifford operator with another Clifford operator."""
@inline function apply_right_row_kernel!(new_lrow, l_tab, row, l_tab_copy, r_tab; phases=true)
    phases && (new_lrow.phase[] = r_tab.phases[row])
    n = nqubits(l_tab)
    for qubit in 1:n
        x,z = r_tab[row,qubit]
        if phases&&x&&z
            new_lrow.phase[] -= 0x1
        end
        if x
            mul_left!(new_lrow, l_tab_copy, qubit, phases=Val(phases))
        end
        if z
            mul_left!(new_lrow, l_tab_copy, qubit+n, phases=Val(phases))
        end
    end
    l_tab[row] = new_lrow
    new_lrow
end


##############################
# Single-qubit gates
##############################

function apply_right!(l::CliffordOperator, r::sHadamard)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sHadamardXY)
    mul_right_ignore_anticomm!(tab(l), r.q, nqubits(l)+r.q)
    apply_right!(l, sY(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sHadamardYZ)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q, r.q)
    apply_right!(l, sZ(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sPhase)
    apply_right!(l, sInvPhase(r.q))
    apply_right!(l, sZ(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvPhase)
    mul_right_ignore_anticomm!(tab(l), r.q, nqubits(l)+r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sX)
    phases(l)[nqubits(l)+r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sY)
    phases(l)[r.q] ⊻= 0x02
    phases(l)[nqubits(l)+r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sZ)
    phases(l)[r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTX)
    apply_right!(l, sInvSQRTX(r.q))
    apply_right!(l, sX(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTX)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q, r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTY)
    phases(l)[nqubits(l)+r.q] ⊻= 0x02
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTY)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    phases(l)[nqubits(l)+r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sCXYZ)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    mul_right_ignore_anticomm!(tab(l), r.q, nqubits(l)+r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sCZYX)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q, r.q)
    apply_right!(l, sX(r.q))
    return l
end

function apply_right!(l::CliffordOperator, ::sId1)
    return l
end


##############################
# Two-qubit gates
##############################

function apply_right!(l::CliffordOperator, r::sSWAP)
    rowswap!(tab(l), nqubits(l)+r.q1, nqubits(l)+r.q2)
    rowswap!(tab(l), r.q1, r.q2)
    return l
end

function apply_right!(l::CliffordOperator, r::sSWAPCX)
    apply_right!(l, sSWAP(r.q1, r.q2))
    apply_right!(l, sCNOT(r.q2, r.q1))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSWAPCX)
    apply_right!(l, sCNOT(r.q2, r.q1))
    apply_right!(l, sSWAP(r.q1, r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sISWAP)
    apply_right!(l, sSWAP(r.q1, r.q2))
    apply_right!(l, sZCZ(r.q1, r.q2))
    apply_right!(l, sPhase(r.q1))
    apply_right!(l, sPhase(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvISWAP)
    apply_right!(l, sSWAP(r.q1, r.q2))
    apply_right!(l, sZCZ(r.q1, r.q2))
    apply_right!(l, sInvPhase(r.q1))
    apply_right!(l, sInvPhase(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sCZSWAP)
    apply_right!(l, sZCZ(r.q2, r.q1))
    apply_right!(l, sSWAP(r.q1, r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sCXSWAP)
    apply_right!(l, sCNOT(r.q2, r.q1))
    apply_right!(l, sSWAP(r.q1, r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sCNOT)
    return apply_right!(l, sZCX(r.q1, r.q2))
end

function apply_right!(l::CliffordOperator, r::sCPHASE)
    apply_right!(l, sZCZ(r.q1, r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sZCX)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q2, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), r.q1, r.q2)
    return l
end

function apply_right!(l::CliffordOperator, r::sZCY)
    apply_right!(l, sHadamardYZ(r.q2))
    apply_right!(l, sZCZ(r.q1, r.q2))
    apply_right!(l, sHadamardYZ(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sZCZ)
    mul_right_ignore_anticomm!(tab(l), r.q2, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), r.q1, nqubits(l)+r.q2)
    return l
end

function apply_right!(l::CliffordOperator, r::sXCX)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q2, r.q1)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q1, r.q2)
    return l
end

function apply_right!(l::CliffordOperator, r::sXCY)
    apply_right!(l, sHadamardXY(r.q2))
    apply_right!(l, sXCX(r.q1, r.q2))
    apply_right!(l, sHadamardXY(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sXCZ)
    return apply_right!(l, sZCX(r.q2, r.q1))
end

function apply_right!(l::CliffordOperator, r::sYCX)
    return apply_right!(l, sXCY(r.q2, r.q1))
end

function apply_right!(l::CliffordOperator, r::sYCY)
    apply_right!(l, sHadamardYZ(r.q1))
    apply_right!(l, sHadamardYZ(r.q2))
    apply_right!(l, sZCZ(r.q1, r.q2))
    apply_right!(l, sHadamardYZ(r.q2))
    apply_right!(l, sHadamardYZ(r.q1))
    return l
end

function apply_right!(l::CliffordOperator, r::sYCZ)
    return apply_right!(l, sZCY(r.q2, r.q1))
end

function apply_right!(l::CliffordOperator, r::sZCrY)
    mul_right_ignore_anticomm!(tab(l), r.q1, r.q2)        
    mul_right_ignore_anticomm!(tab(l), r.q2, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q2, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), r.q1, nqubits(l)+r.q2)
    return l
end

function apply_right!(l::CliffordOperator, r::sInvZCrY)
    mul_right_ignore_anticomm!(tab(l), r.q1, r.q2)
    mul_right_ignore_anticomm!(tab(l), r.q2, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q2, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), r.q1, nqubits(l)+r.q2)
    phases(l)[r.q1] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTZZ)
    apply_right!(l, sInvSQRTZZ(r.q1, r.q2))
    apply_right!(l, sZ(r.q1))
    apply_right!(l, sZ(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTZZ)
    mul_right_ignore_anticomm!(tab(l), r.q1, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), r.q1, nqubits(l)+r.q2)
    mul_right_ignore_anticomm!(tab(l), r.q2, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), r.q2, nqubits(l)+r.q2)
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTXX)
    apply_right!(l, sInvSQRTXX(r.q1, r.q2))
    apply_right!(l, sX(r.q1))
    apply_right!(l, sX(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTXX)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q1, r.q1)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q1, r.q2)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q2, r.q1)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q2, r.q2)    
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTYY)
    apply_right!(l, sInvSQRTYY(r.q1, r.q2))
    apply_right!(l, sY(r.q1))
    apply_right!(l, sY(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTYY)
    mul_right_ignore_anticomm!(tab(l), r.q1, nqubits(l)+r.q1)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q1, nqubits(l)+r.q2)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q1, r.q2)
    mul_right_ignore_anticomm!(tab(l), r.q2, r.q1)
    mul_right_ignore_anticomm!(tab(l), nqubits(l)+r.q2, r.q1)
    mul_right_ignore_anticomm!(tab(l), r.q1, nqubits(l)+r.q1)    
    rowswap!(tab(l), r.q1, nqubits(l)+r.q1)
    rowswap!(tab(l), r.q2, nqubits(l)+r.q2)
    apply_right!(l, sZ(r.q2))
    return l
end


"""Multiply Pauli operators `l * r`, ignoring anticommutation phases (keeping only ±1, not ±i).
See also: [`mul_right!`](@ref)."""
@inline function mul_right_ignore_anticomm!(s::Tableau, m, t::Tableau, i; phases::Val{B}=Val(true)) where B
    extra_phase = mul_right!((@view s.xzs[:,m]), (@view t.xzs[:,i]); phases=phases)
    B && (s.phases[m] = (extra_phase+s.phases[m]+s.phases[i])&0x2)
    s
end
@inline function mul_right_ignore_anticomm!(s::Tableau, m, i; phases::Val{B}=Val(true)) where B
    extra_phase = mul_right!((@view s.xzs[:,m]), (@view s.xzs[:,i]); phases=phases)
    B && (s.phases[m] = (extra_phase+s.phases[m]+s.phases[i])&0x2)
    s
end