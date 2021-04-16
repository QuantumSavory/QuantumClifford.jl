#!/usr/licensed/julia/1.5.0/bin/julia

#= Document me =#
module AncillaCircuits
#=Or maybe =#
using Statistics
using Plots
using Formatting
using QuantumClifford
using LaTeXStrings
using LsqFit
using ProgressMeter
using Base.Threads
using CliffordOps

const MAX_QUBITS = 60;
export get_ancilla_record, get_entropies, follow_trajectory,
get_mutual_sets, get_mutual_end, unitary_layer!, project_layer!

function add_ancilla(s)
    xarr = broadcast(p->push!(xbit(p), 0),s)
    zarr = broadcast(p->push!(zbit(p), 0),s)
    xs = hcat(hcat(xarr...), zeros(Bool,length(xarr[1])))
    zs = hcat(hcat(zarr...), zeros(Bool, length(zarr[1])))
    zs[end,end] = true
    return Stabilizer(vcat(s.phases, [0x0]), xs', zs')
end

# compute mutual info with shrinking patches of ancillae,
# removing the earliest ancillae first
function get_mutual(s, L)
    sys_ent = entropy(s, L)
    mutual = zeros(s.nqubits - L)
    for i in 1:(s.nqubits-L)
        anc_ent = entropy(s, [L+i:s.nqubits;])
        sum_ent = entropy(s, [L+1:L+i-1;])
        mutual[i] = sys_ent + anc_ent - sum_ent
    end
    return mutual
end

# unlike get_mutual, remove ancillae from end
# for total pure state, this will satisfy:
# reverse(mutual) + get_mutual = 2 * sys_ent
# see notebook p 106
function get_mutual_end(s, L)
    sys_ent = entropy(s, L)
    mutual = zeros(s.nqubits - L)
    for i in 1:(s.nqubits-L)
        anc_ent = entropy(s, [L+1:s.nqubits-i;])
        sum_ent = entropy(s, [s.nqubits+1-i:s.nqubits;])
        mutual[i] = sys_ent + anc_ent - sum_ent
    end
    return mutual
end

function get_mutual_sets(s, L;num::Int=1)
    sys_ent = entropy(s, L)
    mutual = zeros(s.nqubits-L-num+1)
    for i in 1:(s.nqubits-L-num+1)
        anc_ent = entropy(s, [L+i:L+i+num-1;])
        sum_ent = entropy(s, [1:L;L+i:L+num+i-1;])
        mutual[i] = sys_ent + anc_ent - sum_ent
    end
    return mutual
end

# each time step: apply layer of 2-qubit random Clifford gates
function unitary_layer!(s, L, parity)
    # I could do this in one loop but just to make things simpler...
    if parity==1
        for i in [1:2:L-1;]
            mat = rand_clifford(2)
            apply!(s,mat,[i,i+1])
        end
    else
        for i in [2:2:L;]
            mat = rand_clifford(2)
            apply!(s,mat,[i,(i+1)%L])
        end
    end
    s
end

# layer of coupling ancillae via CNOT gates
function ancilla_layer!(s, prob, L)
    for i in [1:L;]
        r = rand()
        if r < prob # add an ancilla
            # this doesn't work for >=63 qubits, 
            # so in that case, do the roundabout version with bool arrays
            if s.nqubits < MAX_QUBITS
                s = Stabilizer(vcat(s.phases, [0x0]), s.nqubits + 1, 
                    vcat(s.xzs, [UInt64(0) UInt64(2^(s.nqubits))]))
            else
                s = add_ancilla(s)
            end
            apply!(s, CNOT,[i, s.nqubits])
        end
    end
    return s
end

function project_layer!(s, prob, L)
    for i in [1:L;]
        r = rand()
        if r < prob
            # don't care about phases or anything.
            project!(s, single_z(nqubits(s), i),phases=false,
                     keep_result=false)
        end
    end
end

function dephase_layer!(md::MixedDestabilizer, prob, L)
    for i in [1:L;]
        r = rand()
        if r < prob
            # don't care about phases or anything.
	    # there are some things I need to fix here.
            project!(md, single_z(nqubits(s), i),phases=false,
                     keep_result=false)
        end
    end
end

# compute entropy of system, which is entangled with ancillae
function entropy(s, L::Int)
    if s.nqubits == L
        return 0
    else
        return entropy(s, [L+1:s.nqubits;])
    end
end

# compute entropy over some partition of system
function entropy(s, arr::Array)
    _, rank = canonicalize_rref!(copy(s), arr;phases=false)
    return s.nqubits - length(arr) - rank
end

# Uses different conventions from get_entropies
function get_ancilla_record(probs, initial_runs, num_samples, L; time::Int=4, print::Bool=false, check::Bool=false)
    final_states = Array{Stabilizer}(undef,(length(probs), length(initial_runs), num_samples))
    sample_ancillas = zeros(length(probs), length(initial_runs), num_samples,time*L)
    mutual_infos = Array{Array}(undef, (length(probs), length(initial_runs), num_samples))
    for (prob_i, prob) in enumerate(probs) 
        println("Ancilla rate $(prob), system size $L")
        @showprogress for (init_i, init_run) in enumerate(initial_runs)
            @threads for sample_i in 1:num_samples
                s = one(Stabilizer, L)
                for step in 1:init_run*L
                    unitary_layer!(s,L, step%2)
                end
                print && println(s)
                sample_ancillas[prob_i,init_i, sample_i,1] = 0
                for step in 2:(time*L)
                    unitary_layer!(s,L,step%2)
                    s = ancilla_layer!(s, prob, L)
                    sample_ancillas[prob_i, init_i, sample_i, step] = s.nqubits - L
                end
		if check
		   check_allrowscommute(s) ? nothing : throw(AssertionError(
		   "Rows don't commute in $s, $(prob), $(init_run), $(sample_i)"))
		end

                final_states[prob_i,init_i,sample_i] = s
                mutual_infos[prob_i, init_i, sample_i] = get_mutual(s,L)
            end
        end
    end
    return final_states, sample_ancillas, mutual_infos
end

function dephase_record(probs,num_samples, L; time::Int=4, check::Bool=false,print::Bool=false)
    final_states = Array{Stabilizer}(undef,(length(probs), length(initial_runs), num_samples))
    mutual_infos = Array{Array}(undef, (length(probs), length(initial_runs), num_samples))
    final_entropies = zeros(length(probs), length(initial_runs), num_samples)

    @showprogress for (prob_i, prob) in enumerate(probs)
        @threads for sample_i in 1:num_samples	 
            s = one(Stabilizer, L)

            # now run till entropy is saturated
            while entropy(s,L) < L
            	unitary_layer!(s, L, istep%2)
                s = ancilla_layer!(s,prob,L)
            end
            print && println(s)

            # and now do dephasing measurements
            destab = MixedDestabilizer(s)
            
            for step in 1:time*L
                unitary_layer!(s, L, step%2)
                project_layer!(s, prob, L)
            end
                
            if check
                check_allrowscommute(s) ? nothing : throw(AssertionError(
                "Rows don't commute in $s, $(prob), $(init), $(sample_i)"))
            end
            final_entropies[prob_i, init_i, sample_i] = entropy(s, L)
            final_states[prob_i,init_i,sample_i] = s
            mutual_infos[prob_i, init_i, sample_i] = get_mutual(s,L)
        end
    end
    return final_states, final_entropies, mutual_infos
end

# probably makes no difference how I start because I just want to make it a mixed state...but maybe just in case
function follow_trajectory(probs, initial_runs, num_samples, L; 
                          time::Int=4, check::Bool=false,print::Bool=false)
    final_states = Array{Stabilizer}(undef,(length(probs), length(initial_runs), num_samples))
    mutual_infos = Array{Array}(undef, (length(probs), length(initial_runs), num_samples))
    final_entropies = zeros(length(probs), length(initial_runs), num_samples)
    for (prob_i, prob) in enumerate(probs) 
        println("Ancilla rate $(prob), system size $L")
        @showprogress for (init_i, init_run) in enumerate(initial_runs)
            @threads for sample_i in 1:num_samples
                s = one(Stabilizer, L)
                for step in 1:init_run*L
                   unitary_layer!(s,L,step%2)
                end
                print && println(s)

                # now run till entropy is saturated
                istep = 1
                while entropy(s,L) < L
                   unitary_layer!(s, L, istep%2)
                   s = ancilla_layer!(s,prob,L)
                   istep += 1
                end
                print && println(s)
                print && println("Number of steps to saturation $(istep)")
                # and now do projective measurements
                # num_measure = 0
                for step in 1:time*L
                   unitary_layer!(s, L, step%2)
		   project_layer!(s, prob, L)
                   #num_measure += project_layer!(s,prob,L)
                end
                
                # print && println("Average measurement rate $(num_measure * 1.0/(time * L^2))")
                if check
		   check_allrowscommute(s) ? nothing : throw(AssertionError(
		   "Rows don't commute in $s, $(prob), $(init), $(sample_i)"))
		end
                final_entropies[prob_i, init_i, sample_i] = entropy(s, L)
                final_states[prob_i,init_i,sample_i] = s
                mutual_infos[prob_i, init_i, sample_i] = get_mutual(s,L)
            end
        end
    end
    return final_states, final_entropies, mutual_infos
end

function get_entropies(probs, initial_runs, L; num_samples::Int=100, average::Bool=false, time::Int=64, check::Bool=false,print::Bool=false)
    sample_entropies = zeros(length(probs), length(initial_runs), num_samples, time)
    for (prob_i, prob) in enumerate(probs)
        println("Sampling with p=$(prob), $prob_i")
        @showprogress for (init_i, init_run) in enumerate(initial_runs)
            @threads for sample_i in 1:num_samples
                s = one(Stabilizer, L)
                for step in 1:init_run
                    unitary_layer!(s, L, step%2)
                end
                print && println(s)
                sample_entropies[prob_i,init_i,sample_i,1] = L
                for step in 2:(time)
                    unitary_layer!(s, L, step%2)
                    s = ancilla_layer!(s,prob,L)
                    ent = entropy(s, L)
                    sample_entropies[prob_i,init_i,sample_i,step] = -ent + L
                    if ent == L
                        print && println("Breaking at $step, $(s.nqubits)")
                        break
                    end
                end
		if check
		   check_allrowscommute(s) ? nothing : throw(AssertionError(
		   "Rows don't commute in $s, $(prob), $(init), $(sample_i)"))
		end
            end
        end
    end
    if average
        return [dropdims(mean(sample_entropies, dims=(3)), dims=(3)),
            dropdims(std(sample_entropies, dims=(3))/sqrt(num_samples), dims=(3))]
    else
        return sample_entropies
    end
end

end
