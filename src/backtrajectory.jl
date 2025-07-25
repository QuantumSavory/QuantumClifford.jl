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
        elseif op isa sMX
            push!(results, do_MX!(T, op))
        # elseif op isa sMY
            # push!(results, do_MY!(T, op))
        elseif op isa sMZ
            push!(results, do_MZ!(T, op))
        else
            error("Unsupported operation: $(typeof(op))")
        end
    end

    return results
end


# function do_MY!(T, op::sMY)
#     collapse_y!(T, op.qubit)
#     return eval_y_obs(T, q)
# end

# function collapse_y!(T, q::Int)
#     if is_deterministic_y(T, q)
#         return
#     end

#     apply!(T, sHadamardYZ(q); phases=true)
#     collapse_z!(T, q)
#     apply!(T, sHadamardYZ(q); phases=true)
# end

# function eval_y_obs(T, q::Int)
#     result = T[q]
#     log_i = mul_right!(result, T[nqubits(T)+q])
#     log_i += 1
#     @assert log_i & 1 == 0
#     if log_i & 2
#         result.phase[] 
#     end
#     return result
# end

function do_MX!(T, op::sMX)
    collapse_x!(T, op.qubit)
    return phases(tab(T))[op.qubit] == 0x00 ? 1 : -1
end

# function do_MRX!(T, op::sMRX)
#     collapse_x!(T, op.qubit)
#     result = phases(tab(T))[op.qubit] == 0x00 ? 1 : -1
#     # change the signs to zero
#     return result
# end

function collapse_x!(T, q::Int)
    if is_deterministic_x(T, q)
        return
    end

    apply!(T, sHadamard(q); phases=true)
    collapse_z!(T, q)
    apply!(T, sHadamard(q); phases=true)
end

function do_MZ!(T, op::sMZ)
    collapse_z!(T, op.qubit)
    return phases(tab(T))[op.qubit+nqubits(T)] == 0x00 ? 1 : -1
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


# function backtrajectory(circuit0::Vector{AbstractOperation}, state::AbstractStabilizer)
#     # TODO - Figure out if to use Reset or Gates
#     pushfirst!(circuit0, Reset(state, 1:nqubits(state)))
#     return backtrajectory(circuit0, nqubits(state))
# end

function backtrajectory(circuit::Vector{<:AbstractOperation})
    n = 0
    for op in circuit
        if op isa AbstractSingleQubitOperator
            n = max(n, op.q)
        elseif op isa AbstractTwoQubitOperator
            n = max(n, op.q1, op.q2)
        elseif op isa AbstractMeasurement
            n = max(n, op.qubit)
        else
            error("Unsupported operation: $(typeof(op))")
        end
    end

    return backtrajectory(circuit, n)
end
