module QuantumCliffordPyTesseractDecoderExt

using PyTesseractDecoder: tesseract, stim, np, PythonCall
using SparseArrays
using QuantumClifford
using QECCore

import QuantumClifford.ECC
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, batchdecode, parity_checks

# based on https://github.com/quantumlib/tesseract-decoder/issues/153
function _dem_from_check_matrices(
    check_matrix::SparseMatrixCSC{Bool, Int},
    observables_matrix::SparseMatrixCSC{Bool, Int},
    priors::AbstractVector{<:Real},
)
    num_detectors, num_errors = size(check_matrix)
    num_observables, num_errors2 = size(observables_matrix)
    num_errors == num_errors2 || throw(ArgumentError("check and observable matrices must have the same number of columns"))
    length(priors) == num_errors || throw(DimensionMismatch("expected priors of length $num_errors, got $(length(priors))"))

    dem = stim.DetectorErrorModel()

    for error_col in 1:num_errors
        p = Float64(priors[error_col])
        (0.0 <= p <= 1.0) || throw(DomainError(priors[error_col], "priors must be in the range [0, 1]"))

        targets = PythonCall.Py[]

        for nzidx in nzrange(check_matrix, error_col)
            check_matrix.nzval[nzidx] || continue
            det = check_matrix.rowval[nzidx] - 1
            (0 <= det < num_detectors) || continue
            push!(targets, stim.target_relative_detector_id(det))
        end

        for nzidx in nzrange(observables_matrix, error_col)
            observables_matrix.nzval[nzidx] || continue
            obs = observables_matrix.rowval[nzidx] - 1
            (0 <= obs < num_observables) || continue
            push!(targets, stim.target_logical_observable_id(obs))
        end

        isempty(targets) && continue

        dem.append(stim.DemInstruction(type="error", args=[p], targets=targets))
    end

    return dem
end

function _symptom_to_column_map(
    check_matrix::SparseMatrixCSC{Bool, Int},
    observables_matrix::SparseMatrixCSC{Bool, Int},
)
    num_detectors, num_errors = size(check_matrix)
    num_observables, num_errors2 = size(observables_matrix)
    num_errors == num_errors2 || throw(ArgumentError("check and observable matrices must have the same number of columns"))

    symptom_to_column = Dict{Tuple{Tuple{Vararg{Int}}, Tuple{Vararg{Int}}}, Int}()
    dets = Int[]
    obs = Int[]

    for error_col in 1:num_errors
        empty!(dets)
        empty!(obs)

        for nzidx in nzrange(check_matrix, error_col)
            check_matrix.nzval[nzidx] || continue
            det = check_matrix.rowval[nzidx] - 1
            (0 <= det < num_detectors) || continue
            push!(dets, det)
        end

        for nzidx in nzrange(observables_matrix, error_col)
            observables_matrix.nzval[nzidx] || continue
            o = observables_matrix.rowval[nzidx] - 1
            (0 <= o < num_observables) || continue
            push!(obs, o)
        end

        sort!(dets)
        sort!(obs)
        key = (Tuple(dets), Tuple(obs))
        get!(symptom_to_column, key, error_col)
    end

    return symptom_to_column
end

struct TesseractDecoder <: AbstractSyndromeDecoder
    H
    faults_matrix
    n::Int
    s::Int
    k::Int
    decoder
    error_to_column::Vector{Int}
    nerrors::Int
end

function TesseractDecoder(c; det_beam::Integer=50, priors=nothing, errorrate=nothing)
    H = parity_checks(c)
    s, n = size(H)
    fm = ECC.faults_matrix(H)
    nerrors = 2n

    isnothing(errorrate) || 0 <= errorrate <= 1 || throw(DomainError(errorrate, "`errorrate` must be in the range [0, 1]"))

    prior_vec =
        if !isnothing(priors)
            collect(priors)
        elseif !isnothing(errorrate)
            fill(Float64(errorrate), nerrors)
        else
            fill(0.0001, nerrors)
        end
    length(prior_vec) == nerrors || throw(DimensionMismatch("expected priors of length $nerrors, got $(length(prior_vec))"))

    stab = QuantumClifford.stab_to_gf2(H) # [X | Z]
    check_dense = hcat(@view(stab[:, n+1:end]), @view(stab[:, 1:n])) # [Z | X], maps [X_errors; Z_errors] -> syndrome

    check_sparse = sparse(check_dense)
    fm_sparse = sparse(fm)
    dem = _dem_from_check_matrices(check_sparse, fm_sparse, prior_vec)
    config = tesseract.TesseractConfig(dem=dem, det_beam=Int(det_beam))
    decoder = config.compile_decoder()

    symptom_to_column = _symptom_to_column_map(check_sparse, fm_sparse)
    n_decoder_errors = PythonCall.pyconvert(Int, PythonCall.pybuiltins.len(decoder.errors))
    error_to_column = zeros(Int, n_decoder_errors)
    PythonCall.GIL.@lock begin
        for idx0 in 0:(n_decoder_errors - 1)
            sym = decoder.errors[idx0].symptom
            dets = sort(PythonCall.pyconvert(Vector{Int}, sym.detectors))
            obs = sort(PythonCall.pyconvert(Vector{Int}, sym.observables))
            col = get(symptom_to_column, (Tuple(dets), Tuple(obs)), 0)
            error_to_column[idx0 + 1] = col
        end
    end
    any(==(0), error_to_column) && throw(ArgumentError("unable to map some compiled decoder error indices back to original error columns"))

    k = size(fm, 1) ÷ 2
    return TesseractDecoder(H, fm, n, s, k, decoder, error_to_column, nerrors)
end

parity_checks(d::TesseractDecoder) = d.H

function decode(d::TesseractDecoder, syndrome_sample)
    _s = length(syndrome_sample)
    d.s == _s || throw(ArgumentError(lazy"The syndrome given to `decode` has the wrong dimensions. The syndrome length is $(_s) while it should be $(d.s)"))

    guess = falses(d.nerrors)
    PythonCall.GIL.@lock begin
        idxs = PythonCall.pyconvert(Vector{Int}, d.decoder.decode_to_errors(PythonCall.Py(syndrome_sample).to_numpy()))
        for idx0 in idxs
            (0 <= idx0 < length(d.error_to_column)) || continue
            guess[d.error_to_column[idx0 + 1]] = true
        end
    end
    return guess
end

function batchdecode(d::TesseractDecoder, syndrome_samples)
    samples, _s = size(syndrome_samples)
    d.s == _s || throw(ArgumentError(lazy"The syndromes given to `batchdecode` have the wrong dimensions. The syndrome length is $(_s) while it should be $(d.s)"))

    results = falses(samples, d.nerrors)
    PythonCall.GIL.@lock begin
        for (i, syndrome_sample) in enumerate(eachrow(syndrome_samples))
            idxs = PythonCall.pyconvert(Vector{Int}, d.decoder.decode_to_errors(PythonCall.Py(syndrome_sample).to_numpy()))
            for idx0 in idxs
                (0 <= idx0 < length(d.error_to_column)) || continue
                results[i, d.error_to_column[idx0 + 1]] = true
            end
        end
    end
    return results
end

end
