"""
Simulates measurement results of a Clifford circuit acting on an `n`-qubit |0⟩^⊗n state using the stabilizer tableau backtracking method,
as described by Gidney (2021).

This method incrementally folds operations into an identity tableau by prepending inverses of Clifford gates. Pauli-Z measurements are
resolved by transforming their observables to the initial state; deterministic measurements are directly computed from tableau signs,
while random measurements are simplified and simulated with randomized gate insertions.

Reference:
Gidney, C. (2021). Stim: A fast stabilizer circuit simulator. *Quantum*, 5, 497. https://doi.org/10.22331/q-2021-07-06-497
"""
function backtrajectory(circuit::Vector{<:AbstractOperation}, n::Int)
    T = one(CliffordOperator, n)
    results = Int8[]

    for op in circuit
        if op isa AbstractCliffordOperator
            apply_right!(T, op)
        elseif op isa AbstractMeasurement
            push!(results, do_measure!(T, op))
        elseif op isa sMRX
            push!(results, do_MRX!(T, op))
        elseif op isa sMRY
            push!(results, do_MRY!(T, op))
        elseif op isa sMRZ
            push!(results, do_MRZ!(T, op))
        elseif op isa AbstractReset
            do_reset!(T, op)
        else
            error("Unsupported operation: $(typeof(op))")
        end
    end

    return results
end


function do_measure!(T, op::sMX)
    collapse_x!(T, op.qubit)
    return phases(tab(T))[op.qubit] == 0x00 ? 1 : -1
end

function do_MRX!(T, op::sMRX)
    collapse_x!(T, op.qubit)
    result = phases(tab(T))[op.qubit] == 0x00 ? 1 : -1
    phases(tab(T))[op.qubit] = 0x00
    phases(tab(T))[nqubits(T)+op.qubit] = 0x00
    return result
end

function do_reset!(T, op::sRX)
    collapse_x!(T, op.qubit)
    phases(tab(T))[op.qubit] = 0x00
    phases(tab(T))[nqubits(T)+op.qubit] = 0x00
end

function collapse_x!(T, q::Int)
    if is_deterministic_x(T, q)
        return
    end
    
    apply_right!(T, sHadamard(q))
    collapse_z!(T, q)
    apply_right!(T, sHadamard(q))
end

function do_measure!(T, op::sMY)
    collapse_y!(T, op.qubit)
    return eval_y_obs(T, op.qubit).phase[] == 0x00 ? 1 : -1
end

function do_MRY!(T, op::sMRY)
    collapse_y!(T, op.qubit)
    result = eval_y_obs(T, op.qubit).phase[] == 0x00 ? 1 : -1
    if result == -1
        phases(tab(T))[nqubits(T)+op.qubit] ⊻= 0x02
    end
    return result
end

function do_reset!(T, op::sRY)
    collapse_y!(T, op.qubit)
    if eval_y_obs(T, op.qubit).phase[] != 0x00
        phases(tab(T))[nqubits(T)+op.qubit] ⊻= 0x02
    end
end

function collapse_y!(T, q::Int)
    if is_deterministic_y(T, q)
        return
    end

    apply_right!(T, sHadamardYZ(q))
    collapse_z!(T, q)
    apply_right!(T, sHadamardYZ(q))
end

function eval_y_obs(T, q::Int)
    result = T[q]
    @assert result.phase[] & 0x01 == 0
    og_result_sign = result.phase[]
    mul_right!(result, T[nqubits(T)+q]; phases=Val(true))
    log_i = result.phase[] + 1
    @assert log_i & 0x01 == 0
    if log_i & 2 != 0
        og_result_sign ⊻= 0x02
    end
    result.phase[] = og_result_sign
    return result
end

function do_measure!(T, op::sMZ)
    collapse_z!(T, op.qubit)
    return phases(tab(T))[nqubits(T)+op.qubit] == 0x00 ? 1 : -1
end

function do_MRZ!(T, op::sMRZ)
    collapse_z!(T, op.qubit)
    result = phases(tab(T))[nqubits(T)+op.qubit] == 0x00 ? 1 : -1
    phases(tab(T))[op.qubit] = 0x00
    phases(tab(T))[nqubits(T)+op.qubit] = 0x00
    return result
end

function do_reset!(T, op::sRZ)
    collapse_z!(T, op.qubit)
    phases(tab(T))[op.qubit] = 0x00
    phases(tab(T))[nqubits(T)+op.qubit] = 0x00
end

function collapse_z!(T, q::Int)
    if is_deterministic_z(T, q)
        return
    end

    n = nqubits(T)
    t = tab(T)

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
            apply!(T, sCNOT(pivot, k); phases=true)
        end
    end

    # Swap the now-isolated anti-commuting stabilizer generator for one that commutes with the measurement.
    if getzbit(t, n+q, pivot) == 0
        apply!(T, sHadamard(pivot); phases=true)
    else
        apply!(T, sHadamardYZ(pivot); phases=true)
    end

    # Assign a measurement result.
    if rand(Bool) #TODO: maybe support a RNG
        apply!(T, sX(pivot); phases=true)
    end

    return pivot
end

@inline is_deterministic_x(T, q::Int) = all(getxbytes(T, q) .== 0)
@inline is_deterministic_y(T, q::Int) = all(getxbytes(T, q) .== getxbytes(T, nqubits(T)+q))
@inline is_deterministic_z(T, q::Int) = all(getxbytes(T, nqubits(T)+q) .== 0)

@inline getxbytes(T, r) = tab(T).xzs[1:2:end,r]
@inline getzbytes(T, r) = tab(T).xzs[2:2:end,r]


function do_measure!(T, op::PauliMeasurement)
    if all(iszero.(op.pauli.xz))
        return op.pauli.phase[] & 0x02 == 0x00 ? 1 : -1
    end

    h_xz = []
    h_yz = []
    cnot = []
    meas = 0

    for q in nqubits(op.pauli)
        x, z = op.pauli[q]
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
            push!(cnot, (q, meas))
        end
    end
    @assert meas > 0

    for q in h_xz
        apply_right!(T, sHadamard(q))
    end
    for q in h_yz
        apply_right!(T, sHadamardYZ(q))
    end
    for (q1, q2) in cnot
        apply_right!(T, sCNOT(q1, q2))
    end
    result = do_measure!(T, sMZ(meas))
    for (q1, q2) in reverse(cnot)
        apply_right!(T, sCNOT(q1, q2))
    end
    for q in reverse(h_yz)
        apply_right!(T, sHadamardYZ(q))
    end
    for q in reverse(h_xz)
        apply_right!(T, sHadamard(q))
    end

    return result
end

# function do_set!(T, state)
# end

# function do_reset(T, op::Reset)
# end

# function backtrajectory(circuit::Vector{AbstractOperation}, state::AbstractStabilizer)
#     # TODO - Figure out if to use Reset or Gates
#     pushfirst!(circuit, Reset(state, 1:nqubits(state)))
#     return backtrajectory(circuit, nqubits(state))
# end

function backtrajectory(circuit::Vector{<:AbstractOperation})
    n = maximum(Iterators.flatten(affectedqubits.(circuit)))
    return backtrajectory(circuit, n)
end
