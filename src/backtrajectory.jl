# """Splits some combined operations into their constituent parts."""
# function expand_circuit(circuit::Vector{AbstractOperation})
#     final_circuit = []
#     for op in circuit
#         if op isa sMRX
#         elseif op isa sMRY
#         elseif op isa sMRZ
#         else
#             push!(final_circuit, op)
#         end
#     end
#     return final_circuit
# end

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
function backtrajectory(circuit0::Vector{AbstractOperation}, n::Int)
    # n = nqubits(state)
    # circuit = expand_circuit(circuit0)
    circuit = copy(circuit0)
    T = one(CliffordOperator, n)
    results = []
    
    while !isempty(circuit)
        op = popfirst!(circuit)

        if op isa AbstractCliffordOperator
            T = apply!(inv(CliffordOperator(op, n)), T; phases=true)   # VERY INEFFICIENT - need a fast prepend_inv!
        elseif op isa sMZ
            pivot = 0   
            for c in 1:n
                if getxbit(T.tab, op.qubit+n, c) >= 1   # the Z_q column
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
                    # It’s a Y → Apply S⁻¹ then H
                    apply!(T, sSdg(pivot); phases=true)  #TODO: find s gate
                    apply!(T, sHadamard(pivot); phases=true)
                end
                if rand(Bool)
                    apply!(T, sX(pivot); phases=true)
                end
            end
            push!(results, phases(T.tab)[op.qubit+n] == 0x00 ? 1 : -1)
        elseif op isa Reset
            error("Not implemented: TODO")
        else
            error("Unsupported operation: $(typeof(op))")
        end
    end

    println(T)
    return results
end