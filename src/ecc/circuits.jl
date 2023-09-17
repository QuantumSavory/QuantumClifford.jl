"""Encoding physical qubits into a larger logical code.

Based on [gottesman1997stabilizer](@cite) and [cleve1997efficient](@cite).

The initial physical qubits to be encoded have to be at indices `n-k+1:n`.

!!! info "Encoding circuits are not fault-tolerant"
    Encoding circuits are not fault-tolerant, and thus should not be used in practice.
    Instead, you should measure the stabilizers of the code and the logical observables,
    thus projecting into the code space (which can be fault-tolerant).

!!! warning "Implicit permutation of qubits"
    The canonicalization operation performed on the code may permute the qubits.
    You might need to correct other parts of your code to account for this or
    set `undoperm=true` to add the necessary SWAP gates to undo the permutation.
"""
function naive_encoding_circuit(code; undoperm=false, reportperm=false)
    n = code_n(code)
    k = code_k(code)
    r, permx, permz, md = MixedDestabilizer(code, undoperm=false, reportperm=true);
    circ = QuantumClifford.AbstractOperation[]
    X = logicalxview(md)
    Z = logicalzview(md)
    S = stabilizerview(md)
    for i in 1:k
        for t in 1:n-k
            if X[i,t][1] == true
                push!(circ, sCNOT(n-k+i, t))
            end
        end
    end
    for i in 1:r
        push!(circ, sHadamard(i))
        if S[i,i][2] == true
            push!(circ, sPhase(i))
        end
        for t in 1:n
            if i!=t
                xz = S[i,t]
                g = if xz == (true, true)  # Y
                    sZCY
                elseif xz == (true, false) # X
                    sZCX
                elseif xz == (false, true) && !(i<t<n-k+1) # Z
                    sZCZ
                end
                isnothing(g) || push!(circ, g(i,t))
            end
        end
    end
    for i in 1:n-k # Correct for negative phases in the tableau
        if phases(S)[i]!=0
            if i<=r
                push!(circ, sZ(i))
            else
                push!(circ, sX(i))
            end
        end
    end
    if reportperm
        return circ, permx, permz
    else
        return circ
    end
end

"""Generate the non-fault-tolerant stabilizer measurement cicuit for a given code instance or parity check tableau.

Use the `ancillary_index` and `bit_index` arguments to offset where the corresponding part the circuit starts.
"""
function naive_syndrome_circuit end

function naive_syndrome_circuit(code_type::AbstractECC, ancillary_index=1, bit_index=1)
    naive_syndrome_circuit(parity_checks(code_type), ancillary_index, bit_index)
end

"""Circuit that measures the corresponding PauliOperator by using conditional gates into an ancillary
qubit at index `nqubits(p)+ancillary_index` and stores the measurement result into classical bit `bit_index`."""
function naive_ancillary_paulimeasurement(p::PauliOperator, ancillary_index=1, bit_index=1)
    circuit = AbstractOperation[]
    numQubits = nqubits(p)
    for qubit in 1:numQubits
        if p[qubit] == (1,0)
            push!(circuit, sXCX(qubit, numQubits + ancillary_index))
        elseif p[qubit] == (0,1)
            push!(circuit, sCNOT(qubit, numQubits + ancillary_index))
        elseif p[qubit] == (1,1)
            push!(circuit, sYCX(qubit, numQubits + ancillary_index))
        end
    end
    mz = sMRZ(numQubits + ancillary_index, bit_index)
    push!(circuit, mz)

    return circuit
end

function naive_syndrome_circuit(parity_check_tableau, ancillary_index=1, bit_index=1)
    naive_sc = AbstractOperation[]

    for check in parity_check_tableau
        naive_sc = append!(naive_sc,naive_ancillary_paulimeasurement(check, ancillary_index, bit_index))
        ancillary_index +=1
        bit_index +=1
    end

    return naive_sc
end
