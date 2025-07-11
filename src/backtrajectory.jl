"""Splits combined operations into their constituent parts."""
function expand_circuit(circuit::Vector{<:AbstractOperation})
    final_circuit = AbstractOperation[]
    for op in circuit
        if op isa sMRX
            push!(final_circuit, sMX(op.qubit, op.bit))
            push!(final_circuit, Reset(S"X", op.qubit))
        elseif op isa sMRY
            push!(final_circuit, sMY(op.qubit, op.bit))
            push!(final_circuit, Reset(S"Y", op.qubit))
        elseif op isa sMRZ
            push!(final_circuit, sMZ(op.qubit, op.bit))
            push!(final_circuit, Reset(S"Z", op.qubit))
        else
            push!(final_circuit, op)
        end
    end
    return final_circuit
end


"""the `prepend!` function is used to prepend any quantum operation to unitary Clifford operations"""
function prepend! end
"""the prepend_inv! function is used to prepend the inverse of a quantum operation to unitary Clifford operations"""
function prepend_inv! end

# # SLOW versions of prepend! and prepend_inv! for CliffordOperator
# function prepend!(l::CliffordOperator, r::AbstractCliffordOperator; phases=false)
#     apply!(CliffordOperator(r, nqubits(l)), l; phases=phases)
# end
# function prepend_inv!(l::CliffordOperator, r::AbstractCliffordOperator; phases=false)
#     apply!(inv(CliffordOperator(r, nqubits(l))), l; phases=phases)
# end


# FAST versions of prepend! and prepend_inv! for CliffordOperator - TODO
function prepend!(l::CliffordOperator, r; phases)
    @warn "Slow prepend! operation"
    prepend!(l, CliffordOperator(r, nqubits(l)); phases=phases)
end

function prepend!(l::CliffordOperator, r::CliffordOperator; phases=false)
    nqubits(l)==nqubits(r) || throw(DimensionMismatch("The tableau and the Clifford operator need to act on the same number of qubits. Consider specifying an array of indices as a third argument to the `apply!` function to avoid this error."))
    l_tab = tab(l)
    r_tab = tab(r)
    threadlocal = l.buffer
    new_xzs = Vector{typeof(threadlocal)}(undef, length(l_tab))
    @inbounds for row_r in eachindex(r_tab)
        zero!(threadlocal)
        prepend_row_kernel!(threadlocal, row_r, l_tab, r_tab, phases=phases)
        new_xzs[row_r] = copy(threadlocal)
    end
    @inbounds for row_l in eachindex(l_tab)
        l_tab[row_l] = new_xzs[row_l]
    end
    l
end

@inline function prepend_row_kernel!(new_lrow, row, l_tab, r_tab; phases=true)
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


function prepend_inv!(l::CliffordOperator, r; phases=false)
    @warn "Slow prepend_inv! operation"
    prepend!(l, inv(CliffordOperator(r, nqubits(l))); phases=phases)
end
function prepend_inv!(l::CliffordOperator, r::CliffordOperator; phases=false)
    prepend!(l, inv(r); phases=phases)
end

# Symbolic
function prepend!(l::CliffordOperator, r::sX; phases=false)
    tab(l).phases[nqubits(l)+r.q] ⊻= 0x02
    return l
end
function prepend!(l::CliffordOperator, r::sY; phases=false)
    tab(l).phases[r.q] ⊻= 0x02
    tab(l).phases[nqubits(l)+r.q] ⊻= 0x02
    return l
end
function prepend!(l::CliffordOperator, r::sZ; phases=false)
    tab(l).phases[r.q] ⊻= 0x02
    return l
end

function prepend_inv!(l::CliffordOperator, r::sX; phases=false)
    tab(l).phases[nqubits(l)+r.q] ⊻= 0x02
    return l
end
function prepend_inv!(l::CliffordOperator, r::sY; phases=false)
    error("Not implemented: prepend_inv!(l, r::sY)")
    return l
end
function prepend_inv!(l::CliffordOperator, r::sZ; phases=false)
    tab(l).phases[r.q] ⊻= 0x02
    return l
end


"""
Simulates measurement results of a Clifford circuit acting on an `n`-qubit |0⟩^⊗n state using the stabilizer tableau backtracking method,
as described by Gidney (2021).

This method incrementally folds operations into an identity tableau by prepending inverses of Clifford gates. Pauli-Z measurements are
resolved by transforming their observables to the initial state; deterministic measurements are directly computed from tableau signs,
while random measurements are simplified and simulated with randomized gate insertions.

Notes:
- initial state is assumed to be the all-zero state
- supports Pauli-Z measurement only

Reference:
Gidney, C. (2021). Stim: A fast stabilizer circuit simulator. *Quantum*, 5, 497. https://doi.org/10.22331/q-2021-07-06-497
"""
function backtrajectory(circuit0::Vector{<:AbstractOperation}, n::Int)
    circuit = expand_circuit(circuit0)
    T = one(CliffordOperator, n)
    results = Int8[]
    
    for op in circuit
        if op isa AbstractCliffordOperator
            T = prepend_inv!(T, op; phases=true)
        elseif op isa sMZ
            pivot = 0   
            for c in 1:n
                if getxbit(T.tab, op.qubit+n, c) >= 1
                    if pivot == 0
                        pivot = c
                    else
                        apply!(T, sCNOT(pivot, c); phases=true)
                    end
                end
            end
            if pivot != 0
                if getzbit(T.tab, op.qubit+n, pivot) == 0
                    apply!(T, sHadamard(pivot); phases=true)
                else
                    apply!(T, sHadamardYZ(pivot); phases=true)
                end
                if rand(Bool)
                    apply!(T, sX(pivot); phases=true)
                end
            end
            push!(results, phases(T.tab)[op.qubit+n] == 0x00 ? 1 : -1)
        elseif op isa Reset 
            for i in op.indices
                tab(T).phases[i] = 0x00    # TODO: Check op.resetto
            end
        else
            error("Unsupported operation: $(typeof(op))")
        end
    end

    return results, T
end

# function backtrajectory(circuit0::Vector{AbstractOperation}, state::AbstractStabilizer)
#     # TODO - Figure out if to use Reset or Gates
#     pushfirst!(circuit0, Reset(state, 1:nqubits(state)))
#     return backtrajectory(circuit0, nqubits(state))
# end

# function backtrajectory(circuit0::Vector{AbstractOperation})
#     # TODO: Does not work - need consistency between all AbstractOperation subtypes
#     n = maximum(op.q for op in circuit0)
#     return backtrajectory(circuit0, n)
# end