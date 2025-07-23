"""the `apply_right!` function is used to prepend any quantum operation to unitary Clifford operation"""
function apply_right! end

function apply_right!(l::CliffordOperator, r::AbstractCliffordOperator; phases=true)
    @warn "Slow apply_right! operation: $r"
    apply!(CliffordOperator(r, nqubits(l)), l; phases=phases)
end

# helper
function mul_right_log_i!(l::PauliOperator, r::PauliOperator)
    x = mul_right!(l, r; phases=Val(true))
    x.phase[] = x.phase[] & 0x02
    return x
end

# new symbolics
"""A \"symbolic\" single-qubit SQRTZ. See also : [`SingleQubitOperator`](@ref), [`AbstractSymbolicOperator`](@ref)"""
struct sSQRTZ <: AbstractSingleQubitOperator
    q::Int
    sSQRTZ(q) = if q<0 throw(NoZeroQubit) else new(q) end
end
"""A \"symbolic\" single-qubit InvSQRTZ. See also : [`SingleQubitOperator`](@ref), [`AbstractSymbolicOperator`](@ref)"""
struct sInvSQRTZ <: AbstractSingleQubitOperator
    q::Int
    sInvSQRTZ(q) = if q<0 throw(NoZeroQubit) else new(q) end
end

##############################
# Single-qubit gates
##############################

function apply_right!(l::CliffordOperator, r::sHadamard)
    rowswap!(tab(l), r.q, nqubits(l)+r.q)
    return l
end

function apply_right!(l::CliffordOperator, r::sHadamardXY)
    l[r.q] = mul_right_log_i!(l[r.q], l[nqubits(l)+r.q])
    apply_right!(l, sY(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sHadamardYZ)
    l[nqubits(l)+r.q] = mul_right_log_i!(l[nqubits(l)+r.q], l[r.q])
    apply_right!(l, sZ(r.q))
    return l
end

# function apply_right!(l::CliffordOperator, r::sPhase)
#     return l
# end

# function apply_right!(l::CliffordOperator, r::sInvPhase)
#     return l
# end

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
    l[nqubits(l)+r.q] = mul_right_log_i!(l[nqubits(l)+r.q], l[r.q])
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

function apply_right!(l::CliffordOperator, r::sSQRTZ)
    apply_right!(l, sInvSQRTZ(r.q))
    apply_right!(l, sZ(r.q))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTZ)
    l[r.q] = mul_right_log_i!(l[r.q], l[nqubits(l)+r.q])
    return l
end

# function apply_right!(l::CliffordOperator, r::sCXYZ)
#     return l
# end

# function apply_right!(l::CliffordOperator, r::sCZYX)
#     return l
# end

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

# function apply_right!(l::CliffordOperator, r::sSWAPCX)
#     return l
# end

# function apply_right!(l::CliffordOperator, r::sInvSWAPCX)
#     return l
# end

function apply_right!(l::CliffordOperator, r::sISWAP)
    apply_right!(l, sSWAP(r.q1, r.q2))
    apply_right!(l, sZCZ(r.q1, r.q2))
    apply_right!(l, sSQRTZ(r.q1))
    apply_right!(l, sSQRTZ(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvISWAP)
    apply_right!(l, sSWAP(r.q1, r.q2))
    apply_right!(l, sZCZ(r.q1, r.q2))
    apply_right!(l, sInvSQRTZ(r.q1))
    apply_right!(l, sInvSQRTZ(r.q2))
    return l
end

# function apply_right!(l::CliffordOperator, r::sCZSWAP)
#     return l
# end

# function apply_right!(l::CliffordOperator, r::sCXSWAP)
#     return l
# end

function apply_right!(l::CliffordOperator, r::sCNOT)
    return apply_right!(l, sZCX(r.q1, r.q2))
end

# function apply_right!(l::CliffordOperator, r::sCPHASE)
#     return l
# end

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

# function apply_right!(l::CliffordOperator, r::sZCrY)
#     return l
# end

# function apply_right!(l::CliffordOperator, r::sInvZCrY)
#     return l
# end

function apply_right!(l::CliffordOperator, r::sSQRTZZ)
    apply_right!(l, sInvSQRTZZ(r.q1, r.q2))
    apply_right!(l, sZ(r.q1))
    apply_right!(l, sZ(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTZZ)
    l[r.q1] = mul_right_log_i!(l[r.q1], l[nqubits(l)+r.q1])
    l[r.q1] = mul_right_log_i!(l[r.q1], l[nqubits(l)+r.q2])
    l[r.q2] = mul_right_log_i!(l[r.q2], l[nqubits(l)+r.q1])
    l[r.q2] = mul_right_log_i!(l[r.q2], l[nqubits(l)+r.q2])
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTXX)
    apply_right!(l, sInvSQRTXX(r.q1, r.q2))
    apply_right!(l, sX(r.q1))
    apply_right!(l, sX(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTXX)
    l[nqubits(l)+r.q1] = mul_right_log_i!(l[nqubits(l)+r.q1], l[r.q1])
    l[nqubits(l)+r.q1] = mul_right_log_i!(l[nqubits(l)+r.q1], l[r.q2])
    l[nqubits(l)+r.q2] = mul_right_log_i!(l[nqubits(l)+r.q2], l[r.q1])
    l[nqubits(l)+r.q2] = mul_right_log_i!(l[nqubits(l)+r.q2], l[r.q2])
    return l
end

function apply_right!(l::CliffordOperator, r::sSQRTYY)
    apply_right!(l, sInvSQRTYY(r.q1, r.q2))
    apply_right!(l, sY(r.q1))
    apply_right!(l, sY(r.q2))
    return l
end

function apply_right!(l::CliffordOperator, r::sInvSQRTYY)
    l[r.q1] = mul_right_log_i!(l[r.q1], l[nqubits(l)+r.q1])
    l[nqubits(l)+r.q1] = mul_right_log_i!(l[nqubits(l)+r.q1], l[nqubits(l)+r.q2])
    l[nqubits(l)+r.q1] = mul_right_log_i!(l[nqubits(l)+r.q1], l[r.q2])
    l[r.q2] = mul_right_log_i!(l[r.q2], l[r.q1])
    l[nqubits(l)+r.q2] = mul_right_log_i!(l[nqubits(l)+r.q2], l[r.q1])
    l[r.q1] = mul_right_log_i!(l[r.q1], l[nqubits(l)+r.q1])
    rowswap!(tab(l), r.q1, nqubits(l)+r.q1)
    rowswap!(tab(l), r.q2, nqubits(l)+r.q2)
    apply_right!(l, sZ(r.q2))
    return l
end