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

struct BPOTSDecoder <: AbstractSyndromeDecoder
    original_code       
    H::SparseMatrixCSC{Bool,Int}
    faults_matrix::Matrix{Bool}
    n::Int
    s::Int
    k::Int
    cx::Int
    cz::Int
    bpots_x::LDPCDecoders.BPOTSDecoder
    bpots_z::LDPCDecoders.BPOTSDecoder
end

function BPOTSDecoder(c; errorrate=nothing, maxiter=nothing, T=9, C=2.0)
    # Get stabilizer matrices
    Hx_raw = parity_checks_x(c)
    Hz_raw = parity_checks_z(c)
    H_raw = parity_checks(c)

    # Convert to proper matrices
    if H_raw isa Stabilizer
        H_gf2 = stab_to_gf2(H_raw)
        H = sparse(Bool.(H_gf2))
    else
        H = sparse(Bool.(H_raw))
    end
    
    # Convert X and Z matrices
    if Hx_raw isa Stabilizer
        Hx_gf2 = stab_to_gf2(Hx_raw)
        Hz_gf2 = stab_to_gf2(Hz_raw)
        Hx = sparse(Bool.(Hx_gf2))
        Hz = sparse(Bool.(Hz_gf2))
    else
        Hx = sparse(Bool.(Hx_raw))
        Hz = sparse(Bool.(Hz_raw))
    end
    
    # Get dimensions
    s, n = size(H)
    
    # For quantum codes, determine k
    if c isa Toric || c isa Surface
        k = c isa Toric ? 2 : 1
    else
        k = n - s
    end
    
    cx = size(Hx, 1)
    cz = size(Hz, 1)
    
    # Create fault matrix
    fm = BitMatrix(ones(Bool, s, 2*n))
    
    # Create decoders
    errorrate = something(errorrate, 0.0)
    maxiter = something(maxiter, 200)
    bpots_x = LDPCDecoders.BPOTSDecoder(Hx, errorrate, maxiter; T=T, C=C)
    bpots_z = LDPCDecoders.BPOTSDecoder(Hz, errorrate, maxiter; T=T, C=C)

    # Pass the original code object as the first parameter
    return BPOTSDecoder(c, H, fm, n, s, k, cx, cz, bpots_x, bpots_z)
end

function decode(d::BPOTSDecoder, syndrome_sample::AbstractVector{Bool})
    # Validate input size
    length(syndrome_sample) == d.cx + d.cz || 
        throw(DimensionMismatch("Syndrome length ($(length(syndrome_sample))) does not match expected size ($(d.cx + d.cz))"))
    
    # Split syndrome
    row_x = @view syndrome_sample[1:d.cx]
    row_z = @view syndrome_sample[d.cx+1:d.cx+d.cz]
    
    # Decode both parts
    guess_z, conv_z = LDPCDecoders.decode!(d.bpots_x, Vector(row_x))
    guess_x, conv_x = LDPCDecoders.decode!(d.bpots_z, Vector(row_z))
    
    # Return combined X and Z errors
    return vcat(guess_x, guess_z)
end

function parity_checks(d::BPOTSDecoder)
    return d.H
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
    bfz = LDPCDecoders.BitFlipDecoder(Hz, errorrate, maxiter)

    return BitFlipDecoder(H, fm, n, s, k, cx, cz, bfx, bfz)
end

parity_checks(d::BeliefPropDecoder) = d.H
parity_checks(d::BitFlipDecoder) = d.H

function decode(d::BeliefPropDecoder, syndrome_sample)
    row_x = @view syndrome_sample[1:d.cx]
    row_z = @view syndrome_sample[d.cx+1:d.cx+d.cz]
    guess_z, success = LDPCDecoders.decode!(d.bpdecoderx, row_x)
    guess_x, success = LDPCDecoders.decode!(d.bpdecoderz, row_z)
    return vcat(guess_x, guess_z)
end

function decode(d::BitFlipDecoder, syndrome_sample)
    row_x = @view syndrome_sample[1:d.cx]
    row_z = @view syndrome_sample[d.cx+1:d.cx+d.cz]
    guess_z, success = LDPCDecoders.decode!(d.bfdecoderx, row_x)
    guess_x, success = LDPCDecoders.decode!(d.bfdecoderz, row_z)
    return vcat(guess_x, guess_z)
end

end
