"""the `apply_right!` function is used to right multiply any quantum operation to unitary 
Clifford operation or Pauli product"""
function apply_right! end


# dense_cliffords

function apply_right!(l::CliffordOperator, r::PauliOperator; phases=false)
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The Clifford and Pauli operators need to act on the same number of qubits."))
    tab(l).phases[nqubits(l)+1:end] .⊻= r[1:nqubits(l)]
    tab(l).phases[1:nqubits(l)] .⊻= r[nqubits(l)+1:end]
end

function apply_right!(l::CliffordOperator, r::CliffordOperator; phases=true)
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The two Clifford operators need to act on the same number of qubits."))
    l_tab = tab(l)
    r_tab = tab(r)
    threadlocal = l.buffer
    new_xzs = Vector{typeof(threadlocal)}(undef, length(l_tab))
    @inbounds for row_r in eachindex(r_tab)
        zero!(threadlocal)
        apply_right_row_kernel!(threadlocal, row_r, l_tab, r_tab, phases=phases)
        new_xzs[row_r] = copy(threadlocal)
    end
    @inbounds for row_l in eachindex(l_tab)
        l_tab[row_l] = new_xzs[row_l]
    end
    l
end

@inline function apply_right_row_kernel!(new_lrow, row, l_tab, r_tab; phases=true)
    phases && (new_lrow.phase[] = r_tab.phases[row])
    n = nqubits(l_tab)
    for qubit in 1:n
        x,z = r_tab[row,qubit]
        if phases&&x&&z
            new_lrow.phase[] -= 0x1
        end
        if x
            mul_left!(new_lrow, l_tab, qubit, phases=Val(phases))
        end
        if z
            mul_left!(new_lrow, l_tab, qubit+n, phases=Val(phases))
        end
    end
    new_lrow
end


# symbolic_cliffords

function mul_right_ignore_anticommute!(l::PauliOperator, r::PauliOperator)
    x = mul_right!(l, r; phases=Val(true))
    x.phase[] = x.phase[] & 0x02
    return x
end


##############################
# Single-qubit gates
##############################

function apply_right!(l::CliffordOperator, r::sHadamard)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sHadamardXY)
    l[r.q] = mul_right_ignore_anticommute!(l[r.q], l[nqubits(l)+r.q])
    apply_right!(l, sY(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sHadamardYZ)
    l[nqubits(l)+r.q] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q], l[r.q])
    apply_right!(l, sZ(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sPhase)
    apply_right!(l, sInvPhase(r.q))
    apply_right!(l, sZ(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvPhase)
    l[r.q] = mul_right_ignore_anticommute!(l[r.q], l[nqubits(l)+r.q])
    return l
end

function apply_right!(l::CliffordOperator, r::sX)
    phases(tab(l))[nqubits(l)+r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sY)
    phases(tab(l))[r.q] ⊻= 0x02
    phases(tab(l))[nqubits(l)+r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sZ)
    phases(tab(l))[r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTX)
    apply_right!(l, sInvSQRTX(r.q))
    apply_right!(l, sX(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTX)
    l[nqubits(l)+r.q] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q], l[r.q])
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTY)
    phases(tab(l))[nqubits(l)+r.q] ⊻= 0x02
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTY)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    phases(tab(l))[nqubits(l)+r.q] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sCXYZ)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    l[r.q] = mul_right_ignore_anticommute!(l[r.q], l[nqubits(l)+r.q])
    return l
end

function apply_right!(l::CliffordOperator, r::sCZYX)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    l[nqubits(l)+r.q] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q], l[r.q])
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
    l[nqubits(l)+r.q2] = mul_right!(l[nqubits(l)+r.q2], l[nqubits(l)+r.q1]; phases=Val(true))
    l[r.q1] = mul_right!(l[r.q1], l[r.q2]; phases=Val(true))
    return l
end

function apply_right!(l::CliffordOperator, r::sZCY)
    apply_right!(l, sHadamardYZ(r.q2))
    apply_right!(l, sZCZ(r.q1, r.q2))
    apply_right!(l, sHadamardYZ(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sZCZ)
    l[r.q2] = mul_right!(l[r.q2], l[nqubits(l)+r.q1]; phases=Val(true))
    l[r.q1] = mul_right!(l[r.q1], l[nqubits(l)+r.q2]; phases=Val(true))
    return l
end

function apply_right!(l::CliffordOperator, r::sXCX)
    l[nqubits(l)+r.q2] = mul_right!(l[nqubits(l)+r.q2], l[r.q1]; phases=Val(true))
    l[nqubits(l)+r.q1] = mul_right!(l[nqubits(l)+r.q1], l[r.q2]; phases=Val(true))
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
    l[r.q1] = mul_left!(l[r.q1], l[r.q2])
    l[r.q2] = mul_left!(l[r.q2], l[nqubits(l)+r.q1])
    l[nqubits(l)+r.q2] = mul_left!(l[nqubits(l)+r.q2], l[nqubits(l)+r.q1])
    l[r.q1] = mul_left!(l[r.q1], l[nqubits(l)+r.q2])
    return l
end

function apply_right!(l::CliffordOperator, r::sInvZCrY)
    l[r.q1] = mul_left!(l[r.q1], l[r.q2])
    l[r.q2] = mul_left!(l[r.q2], l[nqubits(l)+r.q1])
    l[nqubits(l)+r.q2] = mul_left!(l[nqubits(l)+r.q2], l[nqubits(l)+r.q1])
    l[r.q1] = mul_left!(l[r.q1], l[nqubits(l)+r.q2])
    phases(tab(l))[r.q1] ⊻= 0x02
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTZZ)
    apply_right!(l, sInvSQRTZZ(r.q1, r.q2))
    apply_right!(l, sZ(r.q1))
    apply_right!(l, sZ(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTZZ)
    l[r.q1] = mul_right_ignore_anticommute!(l[r.q1], l[nqubits(l)+r.q1])
    l[r.q1] = mul_right_ignore_anticommute!(l[r.q1], l[nqubits(l)+r.q2])
    l[r.q2] = mul_right_ignore_anticommute!(l[r.q2], l[nqubits(l)+r.q1])
    l[r.q2] = mul_right_ignore_anticommute!(l[r.q2], l[nqubits(l)+r.q2])
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTXX)
    apply_right!(l, sInvSQRTXX(r.q1, r.q2))
    apply_right!(l, sX(r.q1))
    apply_right!(l, sX(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTXX)
    l[nqubits(l)+r.q1] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q1], l[r.q1])
    l[nqubits(l)+r.q1] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q1], l[r.q2])
    l[nqubits(l)+r.q2] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q2], l[r.q1])
    l[nqubits(l)+r.q2] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q2], l[r.q2])
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTYY)
    apply_right!(l, sInvSQRTYY(r.q1, r.q2))
    apply_right!(l, sY(r.q1))
    apply_right!(l, sY(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTYY)
    l[r.q1] = mul_right_ignore_anticommute!(l[r.q1], l[nqubits(l)+r.q1])
    l[nqubits(l)+r.q1] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q1], l[nqubits(l)+r.q2])
    l[nqubits(l)+r.q1] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q1], l[r.q2])
    l[r.q2] = mul_right_ignore_anticommute!(l[r.q2], l[r.q1])
    l[nqubits(l)+r.q2] = mul_right_ignore_anticommute!(l[nqubits(l)+r.q2], l[r.q1])
    l[r.q1] = mul_right_ignore_anticommute!(l[r.q1], l[nqubits(l)+r.q1])
    rowswap!(tab(l), r.q1, nqubits(l)+r.q1)
    rowswap!(tab(l), r.q2, nqubits(l)+r.q2)
    apply_right!(l, sZ(r.q2))
    return l
end