"""
A brute-force, dynamical-programming based approach on
finding minimal representation (decomposition) of any group elements in terms of generators provided.

We don't need to canonicalize anything because it's single qubit operator and representation is unique.
"""
function decompose_into_generators(e::AbstractSingleQubitOperator, generators::Vector{<:AbstractSingleQubitOperator})
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
        table_prime[SingleQubitOperator(k)] = reverse(v)
    end
    table_prime[SingleQubitOperator(sId1(1))] = [sId1(1)]
    return table_prime
end

const IP_SQRTX_DECOMPOSITION_TABLE = decompose_into_generators(sId1(1), [sInvPhase(1), sSQRTX(1)])

"""Generate a single qubit operator multiplication table"""
function gen_single_qubit_multiplication_table()
    table = Dict{Tuple{SingleQubitOperator, SingleQubitOperator}, SingleQubitOperator}()
    for (a, _) in IP_SQRTX_DECOMPOSITION_TABLE
        for (b, _) in IP_SQRTX_DECOMPOSITION_TABLE
            table[(a, b)] = SingleQubitOperator(a * (b * one(CliffordOperator, 1)))
        end
    end
    return table
end

"""Single qubit operator multiplication table to speed up the `_apply_vop!`"""
const SINGLE_QUBIT_MULTIPLICATION_TABLE = gen_single_qubit_multiplication_table()
