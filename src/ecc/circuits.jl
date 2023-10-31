"""Encoding physical qubits into a larger logical code.

The initial physical qubits to be encoded have to be at indices `n-k+1:n`.

!!! info "Encoding circuits are not fault-tolerant"
    Encoding circuits are not fault-tolerant, and thus should not be used in practice.
    Instead, you should measure the stabilizers of the code and the logical observables,
    thus projecting into the code space (which can be fault-tolerant).

The canonicalization operation performed on the code may permute the qubits (see [canonicalize_gott!](@ref)).
That permutation is corrected for with SWAP gates by default (controlled by the `undoperm` keyword argument).

Based on [gottesman1997stabilizer](@cite) and [cleve1997efficient](@cite),
however it seems the published algorithm has some errors.
Consult the erratum, as well as the more recent [grassl2002algorithmic](@cite) and [grassl2011variations](@cite),
and be aware that this implementation also uses H instead of Z gates.
"""
function naive_encoding_circuit(code; undoperm=true)
    n = code_n(code)
    k = code_k(code)
    md, r, permx, permz = MixedDestabilizer(code, undoperm=false, reportperm=true);
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
    if undoperm
        perm = permx[permz]
        transpositions = perm_to_transpositions(perm)
        for (i,j) in transpositions
            push!(circ, sSWAP(i,j))
        end
    end
    circ
end

function perm_to_transpositions(perm)
    n = length(perm)
    transpositions = Tuple{Int, Int}[]
    for i in n:-1:1
        if perm[i]!=i
            j = findfirst(==(i), perm)
            push!(transpositions, (i, j))
            perm[j] = perm[i]
        end
    end
    return transpositions
end

"""Generate the non-fault-tolerant stabilizer measurement cicuit for a given code instance or parity check tableau.

Use the `ancillary_index` and `bit_index` arguments to offset where the corresponding part the circuit starts.

Returns the circuit, the number of ancillary qubits that were added, and a list of bit indices that will store the measurement results.

See also: [`shor_syndrome_circuit`](@ref)
"""
function naive_syndrome_circuit end

function naive_syndrome_circuit(code_type::AbstractECC, ancillary_index=1, bit_index=1)
    naive_syndrome_circuit(parity_checks(code_type), ancillary_index, bit_index)
end

"""Naive Pauli measurement circuit using a single ancillary qubit.

Not a fault-tolerant circuit, but useful for testing purposes.

Measures the corresponding `PauliOperator` by using conditional gates into an ancillary
qubit at index `nqubits(p)+ancillary_index` and stores the measurement result
into classical bit `bit_index`.

See also: [`naive_syndrome_circuit`](@ref), [`shor_ancillary_paulimeasurement`](@ref)"""
function naive_ancillary_paulimeasurement(p::PauliOperator, ancillary_index=1, bit_index=1)
    circuit = AbstractOperation[]
    n = nqubits(p)
    for qubit in 1:n
        if p[qubit] == (1,0)
            push!(circuit, sXCX(qubit, n + ancillary_index))
        elseif p[qubit] == (0,1)
            push!(circuit, sCNOT(qubit, n + ancillary_index))
        elseif p[qubit] == (1,1)
            push!(circuit, sYCX(qubit, n + ancillary_index))
        end
    end
    p.phase[] == 0 || push!(circuit, sX(n + ancillary_index))
    mz = sMRZ(n + ancillary_index, bit_index)
    push!(circuit, mz)

    return circuit
end

function naive_syndrome_circuit(parity_check_tableau, ancillary_index=1, bit_index=1)
    naive_sc = AbstractOperation[]
    ancillaries = 0
    bits = 0
    for check in parity_check_tableau
        append!(naive_sc,naive_ancillary_paulimeasurement(check, ancillary_index+ancillaries, bit_index+bits))
        ancillaries +=1
        bits +=1
    end

    return naive_sc, ancillaries, bit_index:bit_index+bits-1
end

"""Generate the Shor fault-tolerant stabilizer measurement cicuit for a given code instance or parity check tableau.

Use the `ancillary_index` and `bit_index` arguments to offset where the corresponding part the circuit starts.
Ancillary qubits

Returns:
  - The ancillary cat state preparation circuit.
  - The Shor syndrome measurement circuit.
  - The number of ancillary qubits that were added.
  - The list of bit indices that store the final measurement results.

See also: [`naive_syndrome_circuit`](@ref)
"""
function shor_syndrome_circuit end

function shor_syndrome_circuit(code_type::AbstractECC, ancillary_index=1, bit_index=1)
    shor_syndrome_circuit(parity_checks(code_type), ancillary_index, bit_index)
end

"""Shor's Pauli measurement circuit using a multiple entangled ancillary qubits.

A fault-tolerant circuit.

Measures the corresponding PauliOperator by using conditional gates into multiple ancillary
entangled qubits starting at index `nqubits(p)+ancillary_index`
and stores the measurement result into classical bits starting at `bit_index`.
The final measurement result is the XOR of all the bits.

Returns:
  - The ancillary cat state preparation circuit.
  - The Shor syndrome measurement circuit.
  - One more than the index of the last added ancillary qubit.
  - One more than the index of the last added classical bit.

See also: [`naive_syndrome_circuit`](@ref), [`naive_ancillary_paulimeasurement`](@ref)"""
function shor_ancillary_paulimeasurement(p::PauliOperator, ancillary_index=1, bit_index=1)
    init_ancil_index = ancillary_index
    circuit = AbstractOperation[]
    measurements = AbstractOperation[]
    numQubits = nqubits(p)
    for qubit in 1:numQubits
        if p[qubit] == (1,0)
            push!(circuit, sXCZ(qubit, numQubits + ancillary_index))
        elseif p[qubit] == (0,1)
            push!(circuit, sZCZ(qubit, numQubits + ancillary_index))
        elseif p[qubit] == (1,1)
            push!(circuit, sYCZ(qubit, numQubits + ancillary_index))
        end
        if p[qubit] != (0,0)
            push!(measurements, sMRX(numQubits + ancillary_index, bit_index))
            ancillary_index +=1
            bit_index +=1
        end
    end
    circuit = vcat(circuit,measurements)

    cat_state_circuit = AbstractOperation[]
    push!(cat_state_circuit, sHadamard(numQubits + init_ancil_index))
    for i in (init_ancil_index+1):(ancillary_index -1)
        push!(cat_state_circuit, sCNOT(numQubits + (i - 1), numQubits + i))
    end

    return cat_state_circuit, circuit, ancillary_index, bit_index
end

function shor_syndrome_circuit(parity_check_tableau, ancillary_index=1, bit_index=1)
    shor_sc = AbstractOperation[]
    xor_indices = []
    cat_circuit = AbstractOperation[]
    initial_ancillary_index = ancillary_index

    for check in parity_check_tableau
        cat_circ, circ, new_ancillary_index, new_bit_index = shor_ancillary_paulimeasurement(check, ancillary_index, bit_index)
        push!(xor_indices, collect(bit_index:new_bit_index-1))

        append!(shor_sc, circ)
        append!(cat_circuit, cat_circ)

        ancillary_index = new_ancillary_index
        bit_index = new_bit_index
    end

    final_bits = 0
    for indices in xor_indices
        push!(shor_sc, QuantumClifford.ClassicalXOR(indices, bit_index+final_bits))
        final_bits += 1
    end

    return cat_circuit, shor_sc, ancillary_index-initial_ancillary_index, bit_index:bit_index+final_bits-1
end
