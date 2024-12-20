module QuantumCliffordLDPCDecodersExt

import LDPCDecoders
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, parity_checks

abstract type AbstractLDPCDecoder <: AbstractSyndromeDecoder end

# A common structure to hold shared fields for different decoders
mutable struct GenericLDPCDecoder{D} <: AbstractLDPCDecoder
    H
    faults_matrix
    n
    s
    k
    cx
    cz
    decoderx::D
    decoderz::D
end

# Common constructor for any LDPC decoder (both BitFlip and BeliefProp)
function GenericLDPCDecoder(c, DecoderType; errorrate=nothing, maxiter=nothing)
    Hx = parity_checks_x(c)
    Hz = parity_checks_z(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    cx = size(Hx, 1)
    cz = size(Hz, 1)
    fm = faults_matrix(H)
    isnothing(errorrate) || 0 ≤ errorrate ≤ 1 || error("`errorrate` must be in the range [0, 1].")
    errorrate = isnothing(errorrate) ? 0.0 : errorrate
    maxiter = isnothing(maxiter) ? n : maxiter
    decoderx = DecoderType(Hx, errorrate, maxiter)
    decoderz = DecoderType(Hz, errorrate, maxiter)
    return GenericLDPCDecoder{DecoderType}(H, fm, n, s, k, cx, cz, decoderx, decoderz)
end

function BeliefPropDecoder(c; errorrate=nothing, maxiter=nothing)
    return GenericLDPCDecoder(c, LDPCDecoders.BeliefPropagationDecoder; errorrate, maxiter)
end

function BitFlipDecoder(c; errorrate=nothing, maxiter=nothing)
    return GenericLDPCDecoder(c, LDPCDecoders.BitFlipDecoder; errorrate, maxiter)
end

parity_checks(d::GenericLDPCDecoder) = d.H

function decode(d::GenericLDPCDecoder, syndrome_sample)
    row_x = @view syndrome_sample[1:d.cx]
    row_z = @view syndrome_sample[d.cx+1:d.cx+d.cz]
    guess_z, _ = LDPCDecoders.decode!(d.decoderx, row_x)
    guess_x, _ = LDPCDecoders.decode!(d.decoderz, row_z)
    return vcat(guess_x, guess_z)
end

end
