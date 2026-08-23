module QuantumCliffordPyQDecodersExt

using PyQDecoders: np, sps, ldpc, pm
using PythonCall
using SparseArrays
using QuantumClifford
using QuantumClifford.ECC
using QECCore
import QuantumClifford.ECC: AbstractSyndromeDecoder, decode, batchdecode, parity_checks

abstract type PyBP <: AbstractSyndromeDecoder end

struct PyBeliefPropDecoder <: PyBP # TODO all these decoders have the same fields, maybe we can factor out a common type
    code
    H
    Hx
    Hz
    nx
    nz
    faults_matrix
    pyx
    pyz
end

struct PyBeliefPropOSDecoder <: PyBP # TODO all these decoders have the same fields, maybe we can factor out a common type
    code
    H
    Hx
    Hz
    nx
    nz
    faults_matrix
    pyx
    pyz
end

function PyBeliefPropDecoder(c; maxiter=nothing, bpmethod=nothing, errorrate=nothing)
    Hx = reinterpret(UInt8,collect(parity_matrix_x(c)))
    Hz = reinterpret(UInt8,collect(parity_matrix_z(c)))
    H = parity_checks(c)
    fm = faults_matrix(c)
    max_iter=isnothing(maxiter) ? 0 : maxiter
    bpmethod ∈ (nothing, :productsum, :minsum) || error(lazy"PyBeliefPropDecoder got an unknown belief propagation method argument. `bpmethod` must be one of :productsum, :minsum.")
    bp_method = get(Dict(:productsum => "product_sum", :minsum => "minimum_sum"), bpmethod, "minimum_sum")
    isnothing(errorrate) || 0≤errorrate≤1 || error(lazy"PyBeliefPropDecoder got an invalid error rate argument. `errorrate` must be in the range [0, 1].")
    pyx, pyz = lock(QuantumClifford.ECC._python_decoder_lock) do
        PythonCall.GIL.@lock begin
            error_rate = isnothing(errorrate) ? PythonCall.Py(0.0001) : errorrate
            (
                ldpc.BpDecoder(np.array(Hx); max_iter, bp_method, error_rate),
                ldpc.BpDecoder(np.array(Hz); max_iter, bp_method, error_rate),
            )
        end
    end
    return PyBeliefPropDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

function PyBeliefPropOSDecoder(c; maxiter=nothing, bpmethod=nothing, errorrate=nothing, osdmethod=nothing, osdorder=0)
    Hx = reinterpret(UInt8,collect(parity_matrix_x(c)))
    Hz = reinterpret(UInt8,collect(parity_matrix_z(c)))
    H = parity_checks(c)
    fm = faults_matrix(c)
    max_iter=isnothing(maxiter) ? 0 : maxiter
    bpmethod ∈ (nothing, :productsum, :minsum) || error(lazy"PyBeliefPropOSDecoder got an unknown belief propagation method argument. `bpmethod` must be one of :productsum, :minsum.")
    bp_method = get(Dict(:productsum => "product_sum", :minsum => "minimum_sum"), bpmethod, "minimum_sum")
    isnothing(errorrate) || 0≤errorrate≤1 || error(lazy"PyBeliefPropOSDecoder got an invalid error rate argument. `errorrate` must be in the range [0, 1].")
    isnothing(osdmethod) || osdmethod ∈ (:zeroorder, :exhaustive, :combinationsweep) || error(lazy"PyBeliefPropOSDecoder got an unknown OSD method argument. `osdmethod` must be one of :zeroorder, :exhaustive, :combinationsweep.")
    osd_method = get(Dict(:zeroorder => "OSD_0", :exhaustive => "OSD_E", :combinationsweep => "OSD_CS"), osdmethod, 0)
    osdorder≥0 || error(lazy"PyBeliefPropOSDecoder got an invalid OSD order argument. `osdorder` must be ≥0.")
    osd_order = osdorder
    pyx, pyz = lock(QuantumClifford.ECC._python_decoder_lock) do
        PythonCall.GIL.@lock begin
            error_rate = isnothing(errorrate) ? PythonCall.Py(0.0001) : errorrate
            (
                ldpc.BpOsdDecoder(np.array(Hx); max_iter, bp_method, error_rate, osd_method, osd_order),
                ldpc.BpOsdDecoder(np.array(Hz); max_iter, bp_method, error_rate, osd_method, osd_order),
            )
        end
    end
    return PyBeliefPropOSDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

parity_checks(d::PyBP) = d.H

function decode(d::PyBP, syndrome_sample)
    row_x = @view syndrome_sample[1:d.nx] # TODO figure out a lower-overhead way to move data to python (and make sure the ecc wiki still works after that change -- a lot of the cruft here is because of rare crashes when running the wiki evaluator)
    row_z = @view syndrome_sample[d.nx+1:end]
    return lock(QuantumClifford.ECC._python_decoder_lock) do
        PythonCall.GIL.@lock begin
            guess_z_errors = PythonCall.PyArray(d.pyx.decode(PythonCall.Py(row_x).to_numpy()))
            guess_x_errors = PythonCall.PyArray(d.pyz.decode(PythonCall.Py(row_z).to_numpy()))
            vcat(guess_x_errors, guess_z_errors)
        end
    end
end

struct PyMatchingDecoder <: AbstractSyndromeDecoder # TODO all these decoders have the same fields, maybe we can factor out a common type
    code
    H
    Hx
    Hz
    nx
    nz
    faults_matrix
    pyx
    pyz
end

function PyMatchingDecoder(c; weights=nothing)
    Hx = parity_matrix_x(c) |> collect # TODO keep these sparse
    Hz = parity_matrix_z(c) |> collect
    H = parity_checks(c)
    fm = faults_matrix(c)
    pyx, pyz = lock(QuantumClifford.ECC._python_decoder_lock) do
        PythonCall.GIL.@lock begin
            if isnothing(weights)
                (pm.Matching.from_check_matrix(Hx), pm.Matching.from_check_matrix(Hz))
            else
                (
                    pm.Matching.from_check_matrix(Hx, weights=weights),
                    pm.Matching.from_check_matrix(Hz, weights=weights),
                )
            end
        end
    end
    return PyMatchingDecoder(c, H, Hx, Hz, size(Hx, 1), size(Hz, 1), fm, pyx, pyz)
end

parity_checks(d::PyMatchingDecoder) = d.H

function decode(d::PyMatchingDecoder, syndrome_sample)
    row_x = @view syndrome_sample[1:d.nx]
    row_z = @view syndrome_sample[d.nx+1:end]
    return lock(QuantumClifford.ECC._python_decoder_lock) do
        PythonCall.GIL.@lock begin
            guess_z_errors = PythonCall.PyArray(d.pyx.decode(PythonCall.Py(row_x).to_numpy()))
            guess_x_errors = PythonCall.PyArray(d.pyz.decode(PythonCall.Py(row_z).to_numpy()))
            vcat(guess_x_errors, guess_z_errors)
        end
    end
end

function batchdecode(d::PyMatchingDecoder, syndrome_samples)
    row_x = @view syndrome_samples[:,1:d.nx]
    row_z = @view syndrome_samples[:,d.nx+1:end]
    return lock(QuantumClifford.ECC._python_decoder_lock) do
        PythonCall.GIL.@lock begin
            guess_z_errors = PythonCall.PyArray(d.pyx.decode_batch(row_x))
            guess_x_errors = PythonCall.PyArray(d.pyz.decode_batch(row_z))
            hcat(guess_x_errors, guess_z_errors)
        end
    end
end

end
