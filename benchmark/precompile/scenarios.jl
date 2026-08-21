using Random
using QuantumClifford.ECC: CommutationCheckECCSetup, Shor9, Steane7, TableDecoder,
    code_s, evaluate_decoder, naive_encoding_circuit, naive_syndrome_circuit

function pauli()
    product = P"X" * P"Z"
    tensor = P"X" ⊗ P"Z"
    check(product == P"-iY", "Pauli multiplication failed")
    check(tensor == P"XZ", "Pauli tensor product failed")
    return product, tensor
end

function tableau()
    pauli()
    bell_state = S"-XX
                   +ZZ"
    expected = S"-X_
                 +_Z"
    check(tCNOT * bell_state == expected, "dense Clifford multiplication failed")
    check(apply!(copy(bell_state), tCNOT) == expected, "in-place Clifford application failed")

    state = S"-XXX
              +ZZI
              -IZZ"
    canonicalize!(state)
    mixed = MixedDestabilizer(state)
    projected, anticommuting_index, outcome = project!(mixed, P"ZII")
    check(anticommuting_index != 0, "projection found no anticommuting generator")
    check(isnothing(outcome), "noncommuting projection had a deterministic outcome")
    check(nqubits(projected) == 3, "projection changed the qubit count")
    _, repeated_index, repeated_outcome = project!(copy(projected), P"ZII")
    check(repeated_index == 0, "repeated projection found an anticommuting generator")
    check(repeated_outcome == 0x00, "repeated projection returned the wrong outcome")
    return projected
end

function ecc()
    Random.seed!(0x5143)
    code = Steane7()
    encoding_circuit = naive_encoding_circuit(code)
    syndrome_circuit, ancillaries, syndrome_bits = naive_syndrome_circuit(code)
    frames = pftrajectories(
        [encoding_circuit..., syndrome_circuit...]; trajectories=4, threads=false
    )
    syndromes = measurements(frames)
    check(ancillaries == code_s(code) == 6, "Steane circuit returned the wrong ancillas")
    check(syndrome_bits == 1:6, "Steane circuit returned the wrong syndrome-bit range")
    check(size(syndromes) == (4, 6), "Steane syndrome simulation returned an unexpected shape")
    check(!any(syndromes), "noiseless Steane simulation returned a nonzero syndrome")

    logical_error_rates = evaluate_decoder(
        TableDecoder(Shor9()), CommutationCheckECCSetup(0.001), 32
    )
    check(
        all(rate -> 0.0 <= rate <= 1.0, logical_error_rates),
        "Shor decoder returned an invalid logical error rate",
    )
    return syndromes, logical_error_rates
end

const PRECOMPILE_BENCHMARKS = (
    pauli=pauli,
    tableau=tableau,
    ecc=ecc,
)
