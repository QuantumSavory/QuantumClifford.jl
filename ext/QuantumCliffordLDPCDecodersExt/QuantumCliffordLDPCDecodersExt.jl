module QuantumCliffordLDPCDecodersExt

import LDPCDecoders
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, parity_checks

struct BeliefPropDecoder <: AbstractSyndromeDecoder # TODO all these decoders have the same fields, maybe we can factor out a common type
    H
    faults_matrix
    n
    s
    k
    cx
    cz
    bpdecoderx
    bpdecoderz
end

struct BitFlipDecoder <: AbstractSyndromeDecoder # TODO all these decoders have the same fields, maybe we can factor out a common type
    H
    faults_matrix
    n
    s
    k
    cx
    cz
    bfdecoderx
    bfdecoderz
end

function BeliefPropDecoder(c; errorrate=nothing, maxiter=nothing)
    Hx = parity_checks_x(c)
    Hz = parity_checks_z(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    cx = size(Hx, 1)
    cz = size(Hz, 1)
    fm = faults_matrix(H)

    isnothing(errorrate) || 0≤errorrate≤1 || error(lazy"BeliefPropDecoder got an invalid error rate argument. `errorrate` must be in the range [0, 1].")
    errorrate = isnothing(errorrate) ? 0.0 : errorrate
    maxiter = isnothing(maxiter) ? n : maxiter
    bpx = LDPCDecoders.BeliefPropagationDecoder(Hx, errorrate, maxiter)
    bpz = LDPCDecoders.BeliefPropagationDecoder(Hz, errorrate, maxiter)

    return BeliefPropDecoder(H, fm, n, s, k, cx, cz, bpx, bpz)
end

function BitFlipDecoder(c; errorrate=nothing, maxiter=nothing)
    Hx = parity_checks_x(c)
    Hz = parity_checks_z(c)
    H = parity_checks(c)
    s, n = size(H)
    _, _, r = canonicalize!(Base.copy(H), ranks=true)
    k = n - r
    cx = size(Hx, 1)
    cz = size(Hz, 1)
    fm = faults_matrix(H)

    isnothing(errorrate) || 0≤errorrate≤1 || error(lazy"BitFlipDecoder got an invalid error rate argument. `errorrate` must be in the range [0, 1].")
    errorrate = isnothing(errorrate) ? 0.0 : errorrate
    maxiter = isnothing(maxiter) ? n : maxiter
    bfx = LDPCDecoders.BitFlipDecoder(Hx, errorrate, maxiter)
    bfz = LDPCDecoders.BitFlipDecoder(Hz, errorrate, maxiter)

    return BitFlipDecoder(H, fm, n, s, k, cx, cz, bfx, bfz)
end

parity_checks(d::BeliefPropDecoder) = d.H
parity_checks(d::BitFlipDecoder) = d.H

function decode(d::BeliefPropDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.cx]
    row_z = syndrome_sample[d.cx+1:d.cx+d.cz]
    guess_z, success = LDPCDecoders.decode!(d.bpdecoderx, row_x)
    guess_x, success = LDPCDecoders.decode!(d.bpdecoderz, row_z)
    return vcat(guess_x, guess_z)
end

function decode(d::BitFlipDecoder, syndrome_sample)
    row_x = syndrome_sample[1:d.cx]
    row_z = syndrome_sample[d.cx+1:d.cx+d.cz]
    guess_z, success = LDPCDecoders.decode!(d.bfdecoderx, row_x)
    guess_x, success = LDPCDecoders.decode!(d.bfdecoderz, row_z)
    return vcat(guess_x, guess_z)
end

end
