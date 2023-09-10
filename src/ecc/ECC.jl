module ECC

using QuantumClifford
using QuantumClifford: AbstractOperation
using LinearAlgebra
using Nemo
import QuantumClifford: Stabilizer, MixedDestabilizer

abstract type AbstractECC end

"""The encoding circuit of a given code."""
function encoding_circuit end

"""Parity check tableau of a code."""
function parity_checks end

Stabilizer(c::AbstractECC) = parity_checks(c)
MixedDestabilizer(c::AbstractECC) = MixedDestabilizer(Stabilizer(c))

"""The number of physical qubits in a code."""
function code_n end

"""The number of stabilizer checks in a code."""
function code_s(c::AbstractECC)
    s = length(parity_checks(c))
    return s
end

"""The number of logical qubits in a code."""
function code_k(c::AbstractECC)
    k = code_n(c) - code_s(c)
    return k
end

"""The rate of a code."""
function rate(c::AbstractECC)
    rate = code_k(c)//code_s(c)
    return rate
end

"""Number of physical qubits for a given parity check tableau"""
function code_n(parity_check_tableau)
    return size(parity_check_tableau)[2]
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
  - The cat state preparation circuit.
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
  - The cat state preparation circuit.
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

"""The distance of a code."""
function distance end

"""Parity matrix of a code, given as a stabilizer tableau."""
function parity_matrix(c::AbstractECC)
    paritym = stab_to_gf2(parity_checks(c::AbstractECC))
    return paritym
end

"""Logical X operations of a code."""
function logx_ops(c::AbstractECC)
    md = MixedDestabilizer(parity_checks(c))
    logicalxview(md)
end

"""Logical Z operations of a code."""
function logz_ops(c::AbstractECC)
    md = MixedDestabilizer(parity_checks(c))
    logicalzview(md)
end

"""Error-to-logical-observable map (a.k.a. fault matrix) of a code.

For a code with n physical qubits and k logical qubits this function returns
a 2k × 2n binary matrix O such that
`O[i,j]` is true if the logical observable of index `i`
is flipped by the single physical qubit error of index `j`.
Indexing is such that:

- `O[1:k,:]` is the error-to-logical-X map
- `O[k+1:2k,:]` is the error-to-logical-Z map
- `O[:,1:n]` is the X-physical-error-to-logical map
- `O[n+1:2n,:]` is the Z-physical-error-to-logical map

E.g. for `k=1`, `n=10`, then
if `O[2,5]` is true, then the logical Z observable is flipped by a X₅ error;
and if `O[1,12]` is true, then the logical X observable is flipped by a Z₂ error.

Of note is that there is a lot of freedom in choosing the logical operations!
A logical operator multiplied by a stabilizer operator is still a logical operator.
Similarly there is a different fault matrix for each choice of logical operators.
But once the logical operators are picked, the fault matrix is fixed.

Below we show an example that uses the Shor code. While it is not the smallest code,
it is a convenient choice to showcase the importance of the fault matrix when dealing
with degenerate codes where a correction operation and an error do not need to be the same.

First, consider a single-qubit error, potential correction operations, and their effect on the Shor code:

```jldoctest faults_matrix
julia> using QuantumClifford.ECC: faults_matrix, Shor9

julia> state = MixedDestabilizer(Shor9())
𝒟ℯ𝓈𝓉𝒶𝒷━━━━━
+ Z________
+ ___Z_____
+ _X_______
+ __X______
+ ____X____
+ _____X___
+ ______X__
+ _______X_
𝒳ₗ━━━━━━━━━
+ ______XXX
𝒮𝓉𝒶𝒷━━━━━━━
+ XXX___XXX
+ ___XXXXXX
+ ZZ_______
+ Z_Z______
+ ___ZZ____
+ ___Z_Z___
+ ______Z_Z
+ _______ZZ
𝒵ₗ━━━━━━━━━
+ Z__Z____Z

julia> err_Z₁ = single_z(9,1) # the error we will simulate
+ Z________

julia> cor_Z₂ = single_z(9,2) # the correction operation we will perform
+ _Z_______

julia> err_Z₁ * state # observe that one of the syndrome bits is now flipped
𝒟ℯ𝓈𝓉𝒶𝒷━━━━━
+ Z________
+ ___Z_____
+ _X_______
+ __X______
+ ____X____
+ _____X___
+ ______X__
+ _______X_
𝒳ₗ━━━━━━━━━
+ ______XXX
𝒮𝓉𝒶𝒷━━━━━━━
- XXX___XXX
+ ___XXXXXX
+ ZZ_______
+ Z_Z______
+ ___ZZ____
+ ___Z_Z___
+ ______Z_Z
+ _______ZZ
𝒵ₗ━━━━━━━━━
+ Z__Z____Z

julia> cor_Z₂ * err_Z₁ * state # we are back to a good code state
𝒟ℯ𝓈𝓉𝒶𝒷━━━━━
+ Z________
+ ___Z_____
- _X_______
+ __X______
+ ____X____
+ _____X___
+ ______X__
+ _______X_
𝒳ₗ━━━━━━━━━
+ ______XXX
𝒮𝓉𝒶𝒷━━━━━━━
+ XXX___XXX
+ ___XXXXXX
+ ZZ_______
+ Z_Z______
+ ___ZZ____
+ ___Z_Z___
+ ______Z_Z
+ _______ZZ
𝒵ₗ━━━━━━━━━
+ Z__Z____Z

julia> bad_Z₆Z₉ = single_z(9,6) * single_z(9,9) # a different "correction" operation
+ _____Z__Z

julia> bad_Z₆Z₉ * err_Z₁ * state # the syndrome is trivial, but now we have a logical error
𝒟ℯ𝓈𝓉𝒶𝒷━━━━━
+ Z________
+ ___Z_____
+ _X_______
+ __X______
+ ____X____
- _____X___
+ ______X__
+ _______X_
𝒳ₗ━━━━━━━━━
- ______XXX
𝒮𝓉𝒶𝒷━━━━━━━
+ XXX___XXX
+ ___XXXXXX
+ ZZ_______
+ Z_Z______
+ ___ZZ____
+ ___Z_Z___
+ ______Z_Z
+ _______ZZ
𝒵ₗ━━━━━━━━━
+ Z__Z____Z
```

The success of `cor_Z₂` and the failure of `bad_Z₆Z₉` can be immediately seen through the fault matrix,
as the wrong "correction" does not result in the same logical flips ad the error:

```jldoctest faults_matrix
julia> O = faults_matrix(Shor9())
2×18 BitMatrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1
 1  0  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0

julia> O * stab_to_gf2(err_Z₁)
2-element Vector{Int64}:
 0
 0

julia> O * stab_to_gf2(cor_Z₂)
2-element Vector{Int64}:
 0
 0

julia> O * stab_to_gf2(bad_Z₆Z₉)
2-element Vector{Int64}:
 1
 0
```

While its use in this situation is rather contrived, the fault matrix is incredibly useful
when running large scale simulations in which we want a separate fast error sampling process,
(e.g. with Pauli frames) and a syndrome decoding process, without coupling between them.
We just gather all our syndrome measurement **and logical observables** from the Pauli frame simulations,
and then use them with the fault matrix in the syndrome decoding simulation.
"""
function faults_matrix(c::Stabilizer)
    s, n = size(c)
    k = n-s
    O = falses(2k, 2n)
    md = MixedDestabilizer(c)
    logviews = [logicalxview(md); logicalzview(md)]
    errors = [one(Stabilizer,n; basis=:X);one(Stabilizer,n)]
    for i in 1:2k
        O[i, :] = comm(logviews[i]::PauliOperator, errors) # typeassert for JET
    end
    return O
end

function faults_matrix(c::AbstractECC)
    return faults_matrix(parity_checks(c))
end

# TODO implement isdegenerate

include("./bitflipcode.jl")
include("./shorcode.jl")
include("./steanecode.jl")
include("./cleavecode.jl")
include("./fiveonethreecode.jl")
include("./fivetwotwocode.jl")

function canonicalize_cleve97(checks::Stabilizer)
    d, n = size(checks)
    X0 = (checks |> stab_to_gf2)[:,1:n]'
    Z0 = (checks |> stab_to_gf2)[:,n+1:2n]'

    # Now let's work on getting X1 and Z1
    Z2field = Nemo.ResidueRing(ZZ, 2)
    b = rank(Nemo.matrix(Z2field, X0))

    r = d-b
    k = n-d

    # First let's move the 0 columns to the left:
    nullColumns = []
    notNullColumns = []
    for j in 1:d
        if count(X0[:,j])==0 # TODO make sure this is condition -> null generator
            push!(nullColumns, j)
        else
            push!(notNullColumns, j)
        end
    end

    # Reorder the generators/columns so 0 vectors are in front of X0:
    column_order = vcat(nullColumns, notNullColumns)
    X0_5 = X0[:, column_order]
    Z0_5 = Z0[:, column_order]

    # Second let's perform Gaussian elimination to arrive at X1 and Z1
    bank = []
    for j in (r+1):(r+b)
        i = 1
        STOP = false
        while !STOP 
            if X0_5[i, j] != 1
                i+=1 
                if i > n
                    print("ERROR")
                    STOP = true
                end
            else
                STOP = true
            end
        end
        push!(bank, i)
        for a in (r+1):(r+b)
            if a == j
                continue
            end
            if X0_5[i,a] == 1
                X0_5[:,a] = (X0_5[:,j]+X0_5[:,a]).%2
                Z0_5[:,a] = (Z0_5[:,j]+Z0_5[:,a]).%2
            end
        end
    end

    qubit_order = []
    for i in 1:n
        if i ∉ bank
            push!(qubit_order, i)
        end
    end
    append!(qubit_order, bank)
    X1 = X0_5[qubit_order, :]
    Z1 = Z0_5[qubit_order, :]
    println("Qubit order 1", qubit_order)
    checks = checks[:, qubit_order]

    # Now begins the march towards X2 and Z2, starting with B'
    r1 = rank(Nemo.matrix(Z2field, Z1[1:r+k,1:r]))
    r2 = r - r1

    # First let's move the 0 columns of B to the left:
    nullColumns = []
    notNullColumns = []
    for j in 1:r
        if count(Z1[:,j])==0 # TODO make sure this is condition <-> null generator
            push!(nullColumns, j)
        else
            push!(notNullColumns, j)
        end
    end
    for j in r+1:d
        push!(notNullColumns, j)
    end

    # Reorder the generators/columns so 0 vectors are in front of B:
    column_order = vcat(nullColumns, notNullColumns)
    Z1_5 = Z1[:, column_order]

    # Now Gaussian elimination again, this time over B'
    bank = []
    for j in (r2+1):(r2+r1)
        i = 1
        STOP = false
        while !STOP 
            if Z1_5[i, j] != 1
                i+=1 
                if i > n
                    print("ERROR")
                    STOP = true
                end
            else
                STOP = true
            end
        end
        push!(bank, i)
        for a in (r2+1):(r2+r1)
            if a == j
                continue
            end
            if Z1_5[i,a] == 1
                Z1_5[:,a] = (Z1_5[:,j]+Z1_5[:,a]).%2
            end
        end
    end

    qubit_order = []
    for i in 1:n
        if i ∉ bank && i<= k+r
            push!(qubit_order, i)
        end
    end
    for i in k+r+1:n # rows not in B
        push!(bank, i)
    end
    append!(qubit_order, bank)
    Z2 = Z1_5[qubit_order, :]
    X2 = X1[qubit_order, :] # X is unchanged by operations on b, except for reindexing of qubits
    println("Qubit order 2", qubit_order)
    checks = checks[:, qubit_order]

    # Now time for the final steps before arriving at Xstar Zstar
    B1 = Z2[1:k,1+r2:r2+r1]

    XI = zeros(Bool, k , k)
    for i in 1:k
        XI[i,i] = 1
    end
    Xs = (vcat(XI, zeros(Bool,r2, k), B1', zeros(Bool, b, k)))
    Zs = zeros(Bool, n, k)

    println("k, r2, r1, b: ", k, " ", r2, " ", r1, " ", b)
    Xstar = hcat(Xs,X2)
    Zstar = hcat(Zs,Z2)

    #return X0, X0_5, X1, X2, Xstar, Z0, Z0_5, Z1, Z1_5, Z2, Zstar
    return Xstar, Zstar, Stabilizer(X2', Z2')# TODO Return reordered checks or Stab(X2', Z2')?
end

""" The naive implementation of the encoding circuit by arXiv:quant-ph/9607030 """
function naive_encoding_circuit(checks::Stabilizer)
    d, n = size(checks)
    X0 = (checks |> stab_to_gf2)[:,1:n]'

    Z2field = Nemo.ResidueRing(ZZ, 2)
    b = rank(Nemo.matrix(Z2field, X0))

    r = d-b
    k = n-d

    Xstar, Zstar, standard_tab = canonicalize_cleve97(checks)

    naive_ec = AbstractOperation[]
    for i in (r+k+1):(r+k+b)
        push!(naive_ec, sHadamard(i))
    end

    for j in 1:n
        if count(Xstar[:,j]) == 0 # Secondary generators are ignored
            continue
        end
        if Zstar[j,j]==1
            push!(naive_ec, sZ(j))
        end
        for i in 1:n
            if i == j
                continue
            end
            if Xstar[i,j] == 1 && Zstar[i,j] == 0
                push!(naive_ec, sZCX(j,i))
            elseif Xstar[i,j] == 0 && Zstar[i,j] == 1
                push!(naive_ec, sZCZ(j,i))
            elseif Xstar[i,j] == 1 && Zstar[i,j] == 1
                push!(naive_ec, sZCY(j,i))   
            end
        end
    end

    return naive_ec, standard_tab
end

end #module
