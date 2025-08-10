"""
A quantum state representation for efficient simulation using the stabilizer tableau backtracking 
method [gidney-2021](@cite).

This struct stores the inverse of all Clifford operations applied so far, enabling efficient 
simulation by working backwards from measurements to the initial |0⟩^⊗n state. By conjugating the
current-time observable Pₓ by the inverse Clifford operation we get some observable from the
start of time that is equivalent to measuring Pₓ at the current time.
"""
struct BacktrackingRegister <: AbstractQCState
    inv_circuit::CliffordOperator
    bits::Vector{Bool}
end

BacktrackingRegister(r::BacktrackingRegister) = r
BacktrackingRegister(qbits::Int, mbits::Int=0) = BacktrackingRegister(one(CliffordOperator, qbits), falses(mbits))
# BacktrackingRegister(s::AbstractStabilizer, mbits::Int=0) = BacktrackingRegister(..., falses(mbits))
# BacktrackingRegister(r::Register) = BacktrackingRegister(..., r.bits)

Base.copy(r::BacktrackingRegister) = BacktrackingRegister(copy(r.inv_circuit), copy(r.bits))
Base.:(==)(l::BacktrackingRegister,r::BacktrackingRegister) = l.inv_circuit==r.inv_circuit && l.bits==r.bits
Base.hash(r::BacktrackingRegister, h::UInt) = hash(r.inv_circuit, hash(r.bits, h))

nqubits(r::BacktrackingRegister) = nqubits(r.inv_circuit)
bitview(r::BacktrackingRegister) = r.bits
quantumstate(r::BacktrackingRegister) = apply_inv!(one(Stabilizer, nqubits(r)), r.inv_circuit)
tab(r::BacktrackingRegister) = tab(r.inv_circuit)

tensor(regs::BacktrackingRegister...) = BacktrackingRegister(tensor([r.inv_circuit for r in regs]), [bit for r in regs for bit in r.bits])
# tensor(args::AbstractQCState...) = tensor(BacktrackingRegister.(args)...)

function apply!(r::BacktrackingRegister, op, args...; kwargs...)
    apply_right!(r.inv_circuit, op, args...; kwargs...)
    r
end

function apply!(r::BacktrackingRegister, m::sMX)
    _, res = projectXrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end
function apply!(r::BacktrackingRegister, m::sMY)
    _, res = projectYrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end
function apply!(r::BacktrackingRegister, m::sMZ)
    _, res = projectZrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end
function apply!(r::BacktrackingRegister, m::sMRX)
    _, res = projectXrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    phases(tab(r))[m.qubit] = 0x00
    phases(tab(r))[nqubits(r)+m.qubit] = 0x00
    r
end
function apply!(r::BacktrackingRegister, m::sMRY)
    _, res = projectYrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    if iszero(res)
        phases(tab(r))[nqubits(r)+m.qubit] ⊻= 0x02
    end
    r
end
function apply!(r::BacktrackingRegister, m::sMRZ)
    _, res = projectZrand!(r,m.qubit)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    phases(tab(r))[m.qubit] = 0x00
    phases(tab(r))[nqubits(r)+m.qubit] ⊻= 0x02
    r
end
function apply!(r::BacktrackingRegister, m::PauliMeasurement)
    _, res = projectrand!(r,m.pauli)
    m.bit!=0 && (bitview(r)[m.bit] = !iszero(res))
    r
end

function projectXrand!(r::BacktrackingRegister, m)
    if all(getxrow(tab(r), m) .== 0)
        return r, phases(tab(r))[m]
    end

    apply!(r, sHadamard(m))
    collapse_z!(r.inv_circuit, m)
    apply!(r, sHadamard(m))

    r, phases(tab(r))[m]
end

function projectYrand!(r::BacktrackingRegister, m)
    if all(getxrow(tab(r), m) .== getxrow(tab(r), nqubits(r)+m))
        return r, eval_y_obs(r.inv_circuit, m).phase[]
    end

    apply!(r, sHadamardYZ(m))
    collapse_z!(r.inv_circuit, m)
    apply!(r, sHadamardYZ(m))
    
    r, eval_y_obs(r.inv_circuit, m).phase[]
end

function projectZrand!(r::BacktrackingRegister, m)
    if all(getxrow(tab(r), nqubits(r)+m) .== 0)
        return r, phases(tab(r))[nqubits(r)+m]
    end

    collapse_z!(r.inv_circuit, m)
    
    r, phases(tab(r))[nqubits(r)+m]
end

function projectrand!(r::BacktrackingRegister, pauli)
    if all(iszero.(pauli.xz))
        return pauli.phase[] & 0x02
    end

    h_xz = []
    h_yz = []
    cnot = []
    meas = 0

    for q in nqubits(pauli)
        x, z = pauli[q]
        if x
            if z
                push!(h_yz, q)
            else
                push!(h_xz, q)
            end
        end

        if iszero(meas)
            meas = q
        else
            push!(cnot, q)
        end
    end
    @assert meas > 0

    for q in h_xz
        apply!(r, sHadamard(q))
    end
    for q in h_yz
        apply!(r, sHadamardYZ(q))
    end
    for q1 in cnot
        apply!(r, sCNOT(q1, meas))
    end
    _, res = projectZrand!(r, meas)
    for q1 in reverse(cnot)
        apply!(r, sCNOT(q1, meas))
    end
    for q in reverse(h_yz)
        apply!(r, sHadamardYZ(q))
    end
    for q in reverse(h_xz)
        apply!(r, sHadamard(q))
    end

    r, res
end

# function traceout!(r::BacktrackingRegister, arg)
    # TODO
# end

# function applybranches(r::BacktrackingRegister, op)
    # TODO
# end


function eval_y_obs(c::CliffordOperator, q::Int)
    result = c[q]
    @assert result.phase[] & 0x01 == 0
    og_result_sign = result.phase[]
    mul_right!(result, c[nqubits(c)+q]; phases=Val(true))
    log_i = result.phase[] + 1
    @assert log_i & 0x01 == 0
    if log_i & 2 != 0
        og_result_sign ⊻= 0x02
    end
    result.phase[] = og_result_sign
    return result
end

function collapse_z!(c::CliffordOperator, q::Int)
    n = nqubits(c)
    t = tab(c)

    # Search for any stabilizer generator that anti-commutes with the measurement observable.
    pivot = 1
    while pivot <= n && getxbit(t, n+q, pivot) == 0
        pivot += 1
    end
    if pivot >= n+1
        # No anti-commuting stabilizer generator. Measurement is deterministic.
        return -1
    end

    # Perform partial Gaussian elimination over the stabilizer generators that anti-commute with the measurement.
    # Do this by introducing no-effect-because-control-is-zero CNOTs at the beginning of time.
    for k in pivot+1:n
        if getxbit(t, n+q, k) > 0
            apply!(c, sCNOT(pivot, k); phases=true)
        end
    end

    # Swap the now-isolated anti-commuting stabilizer generator for one that commutes with the measurement.
    if getzbit(t, n+q, pivot) == 0
        apply!(c, sHadamard(pivot); phases=true)
    else
        apply!(c, sHadamardYZ(pivot); phases=true)
    end

    # Assign a measurement result.
    if rand(Bool)
        apply!(c, sX(pivot); phases=true)
    end

    return pivot
end

@inline getxrow(s::Tableau, r::Int) = s.xzs[1:2:end, r]
@inline getzrow(s::Tableau, r::Int) = s.xzs[2:2:end, r]