"""A convenience struct to represent the status of a circuit simulated by [`mctrajectories`](@ref)"""
struct CircuitStatus
    status::Int
end

"""Returned by [`applywstatus!`](@ref) if the circuit simulation should continue."""
const continue_stat = CircuitStatus(0)

"""Returned by [`applywstatus!`](@ref) if the circuit reports a success and there is no undetected error.

See also: [`VerifyOp`](@ref), [`BellMeasurement`](@ref)."""
const true_success_stat = CircuitStatus(1)

"""Returned by [`applywstatus!`](@ref) if the circuit reports a success, but it is a false positive (i.e., there was an undetected error).

See also: [`VerifyOp`](@ref), [`BellMeasurement`](@ref)."""
const false_success_stat = CircuitStatus(2)

"""Returned by [`applywstatus!`](@ref) if the circuit reports a failure.

See also: [`VerifyOp`](@ref), [`BellMeasurement`](@ref)."""
const failure_stat = CircuitStatus(3)

const registered_statuses = ["continue",
                             "true_success",
                             "false_success",
                             "failure"
                             ]

function Base.show(io::IO, s::CircuitStatus)
    print(io, get(registered_statuses,s.status+1,nothing))
    print(io, ":CircuitStatus($(s.status))")
end

"""Used for [`mctrajectories`](@ref)."""
function applywstatus!(state, op)
    apply!(state,op), continue_stat
end

"""Run a single Monte Carlo sample, starting with (and modifying) `state` by applying the given `circuit`. Uses `apply!` under the hood."""
function mctrajectory!(state,circuit)
    for op in circuit
        state, cont = applywstatus!(state, op)
        if cont!=continue_stat
            return state, cont
        end
    end
    return state, continue_stat
end

function countmap(samples::Vector{CircuitStatus}) # A simpler faster version of StatsBase.countmap
    counts = zeros(length(registered_statuses))
    for s in samples
        counts[s.status+1] += 1
    end
    Dict(CircuitStatus(i-1)=>counts[i] for i in eachindex(counts))
end

"""Run multiple Monte Carlo trajectories and report the aggregate final statuses of each.

See also: [`pftrajectories`](@ref), `petrajectories`"""
function mctrajectories(initialstate,circuit;trajectories=500)
    counts = countmap([mctrajectory!(copy(initialstate),circuit)[2] for i in 1:trajectories]) # TODO use threads or at least a generator, but without breaking Polyester
    return counts
end
