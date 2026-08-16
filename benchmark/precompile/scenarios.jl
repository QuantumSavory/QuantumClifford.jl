const total_started_ns = time_ns()
const import_started_ns = total_started_ns
using QuantumClifford
const import_seconds = (time_ns() - import_started_ns) / 1.0e9

using QuantumClifford.ECC: Steane7, naive_encoding_circuit, naive_syndrome_circuit

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
    state = MixedDestabilizer(bell())
    initial = copy(state)
    apply!(state, sHadamard(1))
    apply!(state, sHadamard(1))
    check(state == initial, "two Hadamard gates did not restore the Bell state")

    projected, anticommuting_index, outcome = project!(copy(state), P"Z_")
    check(!isnothing(anticommuting_index), "Bell-state projection found no anticommuting generator")
    check(isnothing(outcome), "Bell-state projection unexpectedly had a deterministic outcome")
    check(nqubits(projected) == 2, "Bell-state projection changed the qubit count")
    return projected
end

function ecc()
    code = Steane7()
    encoding_circuit = naive_encoding_circuit(code)
    syndrome_circuit, _ = naive_syndrome_circuit(code)
    frames = pftrajectories(
        [encoding_circuit..., syndrome_circuit...]; trajectories=4, threads=false
    )
    syndromes = measurements(frames)
    check(size(syndromes) == (4, 6), "Steane syndrome simulation returned an unexpected shape")
    check(!any(syndromes), "noiseless Steane simulation returned a nonzero syndrome")
    return syndromes
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
