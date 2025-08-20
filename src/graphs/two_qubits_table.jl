"""
The subgroup of single qubit Clifford group that commutes with CPHASE.

Specifically, this means
Uᵢ * CPHASE₁₂ = CPHASE₁₂ * Uᵢ for all U ∈ Z_COMMUTATION_SUBGROUP and i ∈ {1,2}
"""
const Z_COMMUTATION_SUBGROUP = [CliffordOperator(x, 1) for x in [sId1(1), sZ(1), sPhase(1), sInvPhase(1)]]

"""
Stabilizers needed to convert between equivalent states to satisfy the commutation constraint (see [`apply_cphase_isolated`](@ref))

!!! info "Large subgroup"

    Luckily, to satisfy the commutation constraint we don't need the whole stabilizer subgroup which is quite big, but only a selected few.
"""
const PAULI_STABILIZERS = Dict(
    # For CPHASE |++>
    true => [
        (sId1(1), sId1(1)), (sX(1), sZ(1)), (sZ(1), sX(1)), (sY(1), sY(1)),
        (CliffordOperator(sHadamard) * CliffordOperator(sZ), CliffordOperator(sZ) * CliffordOperator(sHadamard)),
        (CliffordOperator(sZ) * CliffordOperator(sHadamard), CliffordOperator(sHadamard) * CliffordOperator(sZ)),
        # Not needed for satisfying the commutation constraints
        # but good to keep it here in case we need these stabilizers for anything else
        # (CliffordOperator(sHadamard), CliffordOperator(sHadamard))
    ],
    # For |++>
    false => [
        (sId1(1), sId1(1)),
        (sX(1), sId1(1)),
        (sId1(1), sX(1)),
        (sX(1), sX(1)),
        # (CliffordOperator(sHadamard), CliffordOperator(sHadamard))
    ]
)

"""
Generate a 2-qubit graph state from the graph (represented by whether vertices are connected) and two VOPs
"""
function gen_graph_state(connected::Bool, U1::SingleQubitOperator, U2::SingleQubitOperator)
    if connected
        g = GraphState(S"XZ ZX")
    else
        g = GraphState(S"XI IX")
    end
    vops(g)[1], vops(g)[2] = U1, U2
    return g
end

"""
Given a 2-qubit graph state (characterized by VOPs U₁, U₂ and k∈{0,1}), find VOPs U₁', U₂', m∈{0,1} such that

CPHASE₁₂ U₁ ⊗ U₂ (CPHASE₁₂)ᵏ |++⟩ = U₁' ⊗ U₂' (CPHASE₁₂)ᵐ |++⟩

The the following constraint is satisfied:
- if U1 ∈ Z_COMMUTATION_SUBGROUP then U1' ∈ Z_COMMUTATION_SUBGROUP. Same for U2, U2'.

This is done by applying stabilizer on the resulting state to change VOPs without changing the underlying quantum state.
See also [`PAULI_STABILIZERS`](@ref)
"""
function apply_cphase_isolated(connected::Bool, U1::SingleQubitOperator, U2::SingleQubitOperator)
    g_init = gen_graph_state(connected, U1, U2)
    stab = Stabilizer(copy(g_init))
    apply!(stab, sCPHASE(1, 2))
    g = GraphState(stab)

    found = false
    for s in PAULI_STABILIZERS[connected]
        satisfied = true
        # Compute the VOPs, using statbilizer, that give us an equivalent state
        U_prime = [vops(g)[i] * (s[i] * one(CliffordOperator, 1)) for i in 1:2]
        for qubit_idx in 1:2
            # Check if constraint is satisfied
            if (CliffordOperator(vops(g_init)[qubit_idx], 1) in Z_COMMUTATION_SUBGROUP) && !(U_prime[qubit_idx] in Z_COMMUTATION_SUBGROUP)
                satisfied = false
                break
            end
        end
        if satisfied
            vops(g)[1], vops(g)[2] = SingleQubitOperator(U_prime[1]), SingleQubitOperator(U_prime[2])
            # check again it's indeed equivalent
            if canonicalize!(stab) == canonicalize!(Stabilizer(g))
                found = true
                break
            else
                throw("Stabilizer doesn't seem to stabilize the state. This is not supposed to happen.")
            end
        end
    end
    if !found
        throw("We cannot find any equivalent state matching the commutation constraint.")
    end
    return g
end

"""
Generate a lookup table for isolated two qubits CPHASE gate operation

Specifically, this creates a mapping

(k, U₁, U₂) -> (m, U₁', U₂')

where k, m ∈{0,1} such that

CPHASE₁₂ U₁ ⊗ U₂ (CPHASE₁₂)ᵏ |++⟩ = U₁' ⊗ U₂' (CPHASE₁₂)ᵐ |++⟩

and we guarantee if U1 ∈ Z_COMMUTATION_SUBGROUP then U1' ∈ Z_COMMUTATION_SUBGROUP. Same for U2, U2'.

See also [`ISOLATED_CPHASE_TABLE`](@ref) and [`gen_isolated_cphase_table`](@ref)
"""
function gen_isolated_cphase_table()
    table = Dict{Tuple{Bool,SingleQubitOperator,SingleQubitOperator},GraphState}()

    # Enumerate all possible states of the form (U1 ⊗ U2) (CPHASE)ᵏ |++>, k∈{0,1}
    # There are 24 * 24 * 2 possibilities in total.
    for (U1, _) in IP_SQRTX_DECOMPOSITION_TABLE
        for (U2, _) in IP_SQRTX_DECOMPOSITION_TABLE
            table[(true, U1, U2)] = apply_cphase_isolated(true, U1, U2)
            table[(false, U1, U2)] = apply_cphase_isolated(false, U1, U2)
        end
    end
    return table
end

"""
A mapping

(k, U₁, U₂) -> (m, U₁', U₂')

where k, m ∈{0,1} such that

CPHASE₁₂ U₁ ⊗ U₂ (CPHASE₁₂)ᵏ |++⟩ = U₁' ⊗ U₂' (CPHASE₁₂)ᵐ |++⟩

and we guarantee if U1 ∈ Z_COMMUTATION_SUBGROUP then U1' ∈ Z_COMMUTATION_SUBGROUP. Same for U2, U2'.

See also [`apply_cphase_isolated`](@ref) and [`gen_isolated_cphase_table`](@ref)
"""
const ISOLATED_CPHASE_TABLE = gen_isolated_cphase_table()
