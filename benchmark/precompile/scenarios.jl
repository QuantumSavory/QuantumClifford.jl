const total_started_ns = time_ns()
const import_started_ns = total_started_ns
using QuantumClifford
const import_seconds = (time_ns() - import_started_ns) / 1.0e9

using Random
using QuantumClifford.ECC: CommutationCheckECCSetup, Shor9, Steane7, TableDecoder,
    code_s, evaluate_decoder, naive_encoding_circuit, naive_syndrome_circuit

check(condition, message) = condition || error(message)

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

const SCENARIOS = Dict(
    "pauli" => pauli,
    "tableau" => tableau,
    "ecc" => ecc,
)

length(ARGS) == 1 || error("usage: scenarios.jl SCENARIO")
scenario_name = only(ARGS)
scenario = get(SCENARIOS, scenario_name, nothing)
isnothing(scenario) && error("unknown precompile scenario: $(scenario_name)")

trace_mode = get(ENV, "QC_PRECOMPILE_TRACE", "")
first_result = if trace_mode == "compile"
    @timed Base.@trace_compile scenario()
elseif trace_mode == "dispatch"
    @timed Base.@trace_dispatch scenario()
elseif isempty(trace_mode)
    @timed scenario()
else
    error("QC_PRECOMPILE_TRACE must be empty, compile, or dispatch")
end
total_seconds = (time_ns() - total_started_ns) / 1.0e9
warm_result = @timed scenario()

println(join((
    "RESULT",
    scenario_name,
    import_seconds,
    first_result.time,
    first_result.compile_time,
    first_result.recompile_time,
    total_seconds,
    warm_result.time,
    warm_result.compile_time,
    warm_result.recompile_time,
), '\t'))
