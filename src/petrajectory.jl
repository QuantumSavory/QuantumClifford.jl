"""Compute all possible new states after the application of the given operator. Reports the probability of each one of them. Deterministic (as it reports all branches of potentially random processes), part of the Perturbative Expansion interface."""
function applybranches end

"""Compute all possible new states after the application of the given noise model. Reports the probability of each one of them. Deterministic (as it reports all branches of potentially random processes), part of the Noise interface."""
function applynoise_branches end

#TODO is the use of copy here necessary?
#TODO `copy` does not preserve the `fastcolumn` vs `fastrow` effect
applybranches(state, op; max_order=1) = applybranches(operatordeterminism(typeof(op)), state, op; max_order=max_order)

function applybranches(::DeterministicOperatorTrait, state, op; max_order=1)
    [(applywstatus!(copy(state),op)...,1,0)]
end

function applybranches(::NondeterministicOperatorTrait, state, op::T; max_order=1) where {T<:Union{sMZ,sMX,sMY}}
    s = copy(state)

    if T===sMZ
        d, anticom, res = projectZ!(s, op.qubit)
    elseif T===sMX
        d, anticom, res = projectX!(s, op.qubit)
    elseif T===sMY
        d, anticom, res = projectY!(s, op.qubit)
    else
        error("not reachable")
    end

    if isnothing(res)
        tab(stabilizerview(s)).phases[anticom] = 0x0
        s1 = copy(s)
        tab(stabilizerview(s)).phases[anticom] = 0x2
        s2 = s
        return [(s1, continue_stat, .5, 0), (s2, continue_stat, .5, 0)]
    else 
        return [(s, continue_stat, 1, 0)] end
end

function applybranches(::NondeterministicOperatorTrait, state, op; max_order=1)
    throw(ArgumentError(lazy"""
        You are trying to apply a non-deterministic operator $(typeof(op)) in a perturbative expansion, but this particular operator does not have a `applybranches` method defined for it.
        If this is an operator type that you have defined yourself, please either implement an `applybranches` method for it.
        If you know for a fact that this method is deterministic, then you can opt into the default by giving it a `DeterministicOperatorTrait` trait, e.g. `operatordeterminism(::Type{YourType}) = DeterministicOperatorTrait()`.
        If you were not working with custom operators, please report this as a bug in QuantumClifford.jl -- we probably have not gotten around to expanding this functionality to that particular operator.
    """))
end

function petrajectory(state, circuit; branch_weight=1.0, current_order=0, max_order=1)
    if size(circuit)[1] == 0
        return fill(zero(branch_weight), length(registered_statuses)-1)
    end
    next_op = circuit[1]
    rest_of_circuit = circuit[2:end]

    status_probs = fill(zero(branch_weight), length(registered_statuses)-1)

    for (i,(newstate, status, prob, order)) in enumerate(applybranches(state, next_op, max_order=max_order-current_order))
        if status==continue_stat # TODO is the copy below necessary?
            out_probs = petrajectory(copy(newstate), rest_of_circuit,
                branch_weight=branch_weight*prob, current_order=current_order+order, max_order=max_order)
            status_probs .+= out_probs
        else
            status_probs[status.status] += prob*branch_weight
        end
    end

    return status_probs
end

function petrajectory_keep(state, circuit; branch_weight=1.0, current_order=0, max_order=1) # TODO a lot of repetition with petrajectory - dry out
    A = Accumulator{Tuple{typeof(state),CircuitStatus},typeof(branch_weight)}
    dict = A() 
    if size(circuit)[1] == 0
        DataStructures.inc!(dict, (state,continue_stat), branch_weight)
        return dict
    end

    next_op = circuit[1]
    rest_of_circuit = circuit[2:end]

    dict = A()

    for (i,(newstate, status, prob, order)) in enumerate(applybranches(state, next_op, max_order=max_order-current_order))
        if status==continue_stat # TODO is the copy below necessary?
            out_dict = petrajectory_keep(copy(newstate), rest_of_circuit,
                branch_weight=branch_weight*prob, current_order=current_order+order, max_order=max_order)
            DataStructures.merge!(dict, out_dict)
        else
            DataStructures.inc!(dict, (newstate,status), prob*branch_weight)
        end
    end

    return dict
end

"""Run a perturbative expansion to a given order. This is the main public function for the perturbative expansion approach.
Currently only has defined behavior for Clifford operators(AbstractCliffordOperator) and single-qubit X, Y and Z measurements(sMX, sMY, sMZ)
If using the measurements, set keepstates to true. When keepstates is false, will return 0 for true success probability.
See also: [`pftrajectories`](@ref), [`mctrajectories`](@ref)"""
function petrajectories(initialstate, circuit; branch_weight=1.0, max_order=1, keepstates::Bool=false)
    if keepstates
        return petrajectory_keep(initialstate, circuit; branch_weight=branch_weight, current_order=0, max_order=max_order)
    else
        status_probs = petrajectory(initialstate, circuit; branch_weight=branch_weight, current_order=0, max_order=max_order)
        return Dict([CircuitStatus(i)=>status_probs[i] for i in eachindex(status_probs)])
    end
end
