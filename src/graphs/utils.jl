# Code and tables useful to graph state simulation.

# A brute-force, dynamical-programming based approach on
# finding minimal representation (decomposition) of any group elements in terms of generators provided.
#
# We don't need to canonicalize anything because it's single qubit operator and representation is unique.
function decompose(e::AbstractSingleQubitOperator, generators::Vector{<:AbstractSingleQubitOperator})
    e = CliffordOperator(e, 1)

    # Initialize the table with generators. Table shouldn't include identity `e`
    # because composing with identity gives nothing new.
    table = Dict{CliffordOperator,Vector{AbstractSingleQubitOperator}}()
    for g in generators
        push!(table, CliffordOperator(g, 1) => [g])
    end

    updated = true
    while updated
        updated = false
        k = copy(keys(table))
        for a in k
            for b in k
                # compose a new element and see if it's already in the table
                if (a * b) in k
                    if length(table[a]) + length(table[b]) < length(table[a*b])
                        table[a*b] = append!(copy(table[a]), table[b])
                        updated = true
                    end
                else
                    # Don't store identity in the table
                    if a * b != e
                        table[a*b] = append!(copy(table[a]), table[b])
                        updated = true
                    end
                end
            end
        end
    end

    # transform the key to SingleQubitOperator so it's easier to lookup
    table_prime = Dict{SingleQubitOperator,Vector{AbstractSingleQubitOperator}}()
    for (k, v) in table
        table_prime[SingleQubitOperator(k)] = v
    end
    return table_prime
end

# The subgroup of single qubit Clifford group that commutes with CPHASE.
const Z_COMMUTATION_SUBGROUP = [CliffordOperator(x, 1) for x in [sId1(1), sZ(1), sPhase(1), sInvPhase(1)]]

# Stabilizers needed to convert between equivalent states to satisfy the commutation constraint (see cphase_isolated)
# NOTE: Luckily, we don't need the whole stabilizer subgroup which is quite big, but only a selected few.
const PAULI_STABILIZERS = Dict(
    # CPHASE |++>
    true => [
        (sId1(1), sId1(1)), (sX(1), sZ(1)), (sZ(1), sX(1)), (sY(1), sY(1)),
        (CliffordOperator(sHadamard) * CliffordOperator(sZ), CliffordOperator(sZ) * CliffordOperator(sHadamard)),
        (CliffordOperator(sZ) * CliffordOperator(sHadamard), CliffordOperator(sHadamard) * CliffordOperator(sZ)),
        #  (CliffordOperator(sHadamard), CliffordOperator(sHadamard))
    ],
    # |++>
    false => [
        (sId1(1), sId1(1)),
        (sX(1), sId1(1)),
        (sId1(1), sX(1)),
        (sX(1), sX(1)),
        #  (CliffordOperator(sHadamard), CliffordOperator(sHadamard))
    ]
)

function gen_graph_state(connected::Bool, U1::SingleQubitOperator, U2::SingleQubitOperator)
    if connected
        g = GraphState(S"XZ ZX")
    else
        g = GraphState(S"XI IX")
    end
    vops(g)[1], vops(g)[2] = U1, U2
    return g
end

# Convert the graphstate to stabilizer state and perform gate there.
# The resulting state satisfy the following constraint:
# - Let U1, U2 be VOPs of the initial state, U1', U2' be the VOPs of the state after CPHASE
# - We guarantee U1' ∈ Z_COMMUTATION_SUBGROUP if U1 ∈ Z_COMMUTATION_SUBGROUP. Same for U2, U2'.
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

# Generate a lookup table for isolated two qubits CPHASE gate operation
# keys are:
# - whether edge connected
# - first VOP
# - second VOP
# CPHASE acts with control bit on the first qubit.
function gen_isolated_cphase_table()
    table = Dict{Tuple{Bool,SingleQubitOperator,SingleQubitOperator},GraphState}()


    t = copy(IP_SQRTX_DECOMPOSITION_TABLE)
    t[SingleQubitOperator(sId1(1))] = [sId1(1)]

    # Enumerate all possible states of the form (U1 ⊗ U2) (CPHASE)ᵏ |++>, k∈{0,1}
    # There are 24 * 24 * 2 possibilities in total.
    for (U1, _) in t
        for (U2, _) in t
            table[(true, U1, U2)] = apply_cphase_isolated(true, U1, U2)
            table[(false, U1, U2)] = apply_cphase_isolated(false, U1, U2)
        end
    end
    return table
end

# Tables below are precompiled
const IP_SQRTX_DECOMPOSITION_TABLE = decompose(sId1(1), [sInvPhase(1), sSQRTX(1)])
const ISOLATED_CPHASE_TABLE = gen_isolated_cphase_table()
