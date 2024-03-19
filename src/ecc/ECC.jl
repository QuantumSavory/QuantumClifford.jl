module ECC

using LinearAlgebra
using QuantumClifford
using QuantumClifford: AbstractOperation, AbstractStabilizer, Stabilizer
import QuantumClifford: Stabilizer, MixedDestabilizer
using DocStringExtensions
using Combinatorics: combinations
using SparseArrays: sparse
using Statistics: std
using Nemo: ZZ, residue_ring, matrix

abstract type AbstractECC end

export parity_checks, parity_checks_x, parity_checks_z, iscss,
    code_n, code_s, code_k, rate, distance,
    isdegenerate, faults_matrix,
    naive_syndrome_circuit, shor_syndrome_circuit, naive_encoding_circuit,
    RepCode,
    CSS,
    Shor9, Steane7, Cleve8, Perfect5, Bitflip3,
    Toric, Gottesman,
    evaluate_decoder,
    CommutationCheckECCSetup, NaiveSyndromeECCSetup, ShorSyndromeECCSetup,
    TableDecoder,
    BeliefPropDecoder,
    PyBeliefPropDecoder, PyBeliefPropOSDecoder, PyMatchingDecoder

"""Parity check tableau of a code.

See also: [`parity_checks_x`](@ref) and [`parity_checks_z`](@ref)"""
function parity_checks end

"""Parity check boolean matrix of a code (only the X entries in the tableau, i.e. the checks for Z errors).

Only CSS codes have this method.

See also: [`parity_checks`](@ref)"""
function parity_checks_x(code::AbstractECC)
    throw(lazy"Codes of type $(typeof(code)) do not have separate X and Z parity checks, either because they are not a CSS code and thus inherently do not have separate checks, or because its separate checks are not yet implemented in this library.")
end

"""Parity check boolean matrix of a code (only the Z entries in the tableau, i.e. the checks for X errors).

Only CSS codes have this method.

See also: [`parity_checks`](@ref)"""
function parity_checks_z(code::AbstractECC)
    throw(lazy"Codes of type $(typeof(code)) do not have separate X and Z parity checks, either because they are not a CSS code and thus inherently do not have separate checks, or because its separate checks are not yet implemented in this library.")
end

function iscss(::Type{T}) where T<:AbstractECC
    return false
end

function iscss(c::AbstractECC)
    return iscss(typeof(c))
end

parity_checks(s::Stabilizer) = s
Stabilizer(c::AbstractECC) = parity_checks(c)
MixedDestabilizer(c::AbstractECC; kwarg...) = MixedDestabilizer(Stabilizer(c); kwarg...)

"""The number of physical qubits in a code."""
function code_n end

code_n(c::AbstractECC) = code_n(parity_checks(c))

code_n(s::Stabilizer) = nqubits(s)

"""The number of stabilizer checks in a code."""
function code_s end

code_s(c::AbstractECC) = code_s(parity_checks(c))
code_s(s::Stabilizer) = length(s)

"""The number of logical qubits in a code."""
code_k(c) = code_n(c) - code_s(c)

"""The rate of a code."""
function rate(c)
    rate = code_k(c)//code_n(c)
    return rate
end


"""The distance of a code."""
function distance end

"""Parity matrix of a code, given as a stabilizer tableau."""
function parity_matrix(c::AbstractECC)
    paritym = stab_to_gf2(parity_checks(c::AbstractECC))
    return paritym
end

"""Logical X operations of a code."""
function logx_ops(c)
    md = MixedDestabilizer(parity_checks(c))
    logicalxview(md)
end

"""Logical Z operations of a code."""
function logz_ops(c)
    md = MixedDestabilizer(parity_checks(c))
    logicalzview(md)
end

"""Error-to-logical-observable map (a.k.a. fault matrix) of a code.

For a code with n physical qubits and k logical qubits this function returns
a 2k × 2n binary matrix O such that
`O[i,j]` is true if the logical observable of index `i`
is flipped by the single physical qubit error of index `j`.
Indexing is such that:

- `O[1:k,:]` is the error-to-logical-X-observable map (logical X observable, i.e. triggered by logical Z errors)
- `O[k+1:2k,:]` is the error-to-logical-Z-observable map
- `O[:,1:n]` is the X-physical-error-to-logical-observable map
- `O[n+1:2n,:]` is the Z-physical-error-to-logical-observable map

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
    md = MixedDestabilizer(c)
    s, n = size(c)
    r = rank(md)
    k = n - r
    k == n-s || @warn "`faults_matrix` was called on an ECC that has redundant rows (is rank-deficient). `faults_matrix` corrected for that, however this is a frequent source of mistakes and inefficiencies. We advise you remove redundant rows from your ECC."
    O = falses(2k, 2n)
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

include("circuits.jl")
include("decoder_pipeline.jl")

include("codes/util.jl")
include("codes/classical_codes.jl")
include("codes/css.jl")
include("codes/bitflipcode.jl")
include("codes/fivequbit.jl")
include("codes/steanecode.jl")
include("codes/shorcode.jl")
include("codes/clevecode.jl")
include("codes/toric.jl")
include("codes/gottesman.jl")
end #module
