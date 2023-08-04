module ECC

using LinearAlgebra
using QuantumClifford
using QuantumClifford: AbstractOperation, AbstractStabilizer
import QuantumClifford: Stabilizer, MixedDestabilizer
using DocStringExtensions
using Combinatorics: combinations

abstract type AbstractECC end

export Shor9, Steane7, Steane5, Cleve8,
    parity_checks, naive_syndrome_circuit, encoding_circuit,
    code_n, code_s, code_k, rate, distance,
    isdegenerate, faults_matrix

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

"""
$TYPEDSIGNATURES

Check if the code is degenerate with respect to a given set of error or with respect to all
"up to d physical-qubit" errors (defaulting to d=1).

```jldoctest
julia> using QuantumClifford.ECC

julia> isdegenerate(Shor9(), [single_z(9,1), single_z(9,2)])
true

julia> isdegenerate(Shor9(), [single_z(9,1), single_x(9,1)])
false

julia> isdegenerate(Steane7(), 1)
false

julia> isdegenerate(Steane7(), 2)
true
```
"""
function isdegenerate end

isdegenerate(c::AbstractECC, args...) = isdegenerate(parity_checks(c), args...)
isdegenerate(c::AbstractStabilizer, args...) = isdegenerate(stabilizerview(c), args...)

function isdegenerate(H::Stabilizer, errors) # Described in https://quantumcomputing.stackexchange.com/questions/27279
    syndromes = comm.((H,), errors) # TODO This can be optimized by having something that always returns bitvectors
    return length(Set(syndromes)) != length(errors)
end

function isdegenerate(H::Stabilizer, d::Int=1)
    n = nqubits(H)
    errors = [begin p=zero(PauliOperator,n); for i in bits; p[i]=op; end; p end for bits in combinations(1:n, d) for op in ((true,false), (false,true))]
    return isdegenerate(H, errors)
end

"""The standardized logical tableau of a code by [PhysRevA.56.76](@cleve1997efficient)"""
function standard_tab_gott(c::AbstractECC; undoperm= true)
    r, permx, permz, destab = MixedDestabilizer(parity_checks(c); undoperm=undoperm, reportperm=true) # originally undoperm here = false and the comment below is uncommented
    n, k = code_n(c),code_k(c)
    standard_tab = vcat(stabilizerview(destab), logicalxview(destab))
    standard_tab = stab_to_gf2(standard_tab)
    standard_tab = hcat(transpose(reverse(standard_tab[1:n, 1:n])), transpose(reverse(standard_tab[1:n, n+1:2*n])))
    return r, standard_tab
end

function standard_code_tab(c::AbstractECC)
    n, s, k = code_n(c), code_s(c), code_k(c)
    standard_tab = stab_to_gf2(stabilizerview(MixedDestabilizer(parity_checks(c), undoperm=false)))
    md = MixedDestabilizer(parity_checks(c), undoperm=false)
    XZᵗ = stab_to_gf2(stabilizerview(md))
    X₂ = reverse(XZᵗ[1:s,1:n]', dims=(1,2))
    Z₂ = reverse(XZᵗ[1:s,n+1:2n]', dims=(1,2))
    X = falses(n, n)
    Z = falses(n, n)
    X[:, k+1:end] .= X₂
    Z[:, k+1:end] .= Z₂
    X[:, 1:k] .= Xˡ
    return X, Z
    # TODO The permutations need to be reverted at some point, otherwise the generated circuit will have permuted qubits. It is not clear to me whether it is more convenient to revert the permutation here or in naive_encoding_circuit
    # TODO This function does not seem to actually be necessary. It seems like iterating over the `MixedDestabilizer` object is enough. Out of convenience, let's keep if for the moment, but also keep this TODO
end

""" The naive implementation of the encoding circuit by arXiv:quant-ph/9607030 """
function naive_encoding_circuit(c::AbstractECC)
    n= code_n(c)
    r, standard_tab = standard_tab_gott(c)

    naive_ec = AbstractOperation[]
    # Applying the hadamard gate to the last r qubits
    for i in n:-1:n-r+1
        push!(naive_ec, sHadamard(i))
    end
    for i in 1:n
        if standard_tab[i, i] ==1
            for j in 1:n
                if j == i continue end
                gate = (standard_tab[j, i], standard_tab[j, i+n])
                if gate == (1,0)
                    push!(naive_ec, sCNOT(j, i))
                elseif gate == (0,1)
                    push!(naive_ec, sXCZ(j, i))
                elseif gate == (1,1)
                    push!(naive_ec, sXCY(j, i))
                end
            end
        end
    end
    # if undoperm
    #     standard_tab = standard_tab[:,invperm(permx[permz])]
    # end
    return naive_ec
end

include("./bitflipcode.jl")
include("./shorcode.jl")
include("./steanecode.jl")
include("./clevecode.jl")

end #module
