using Graphs: inneighbors, has_edge, rem_edge!, add_edge!

"""
Toggle the edge of a graph between vertices v1 and v2.

Uses `rem_edge!` and `add_edge!`
"""
function toggle_edge!(g, v1::Int, v2::Int)
    if has_edge(g, v1, v2)
        rem_edge!(g, v1, v2)
    else
        add_edge!(g, v1, v2)
    end
    return g
end

"""
Perform Local Complementation repeatedly to transfer the VOP of src to the target.

See also [`local_comp!`](@ref)
"""
function remove_vop!(gs::GraphState, src::Int, target::Int)
    src_seq = IP_SQRTX_DECOMPOSITION_TABLE[vops(gs)[src]]

    if !has_edge(graph(gs), src, target)
        throw("Target qubit must be a neighbor of the source qubit")
    end

    # Perform the local complementation operation to remove generators in the decomposition one-by-one
    for v in src_seq
        if v isa sInvPhase
            # perform local_comp on the target qubit to remove sInvPhase
            local_comp!(gs, target)
        elseif v isa sSQRTX
            # perform local_comp on the source qubit to remove sSQRTX
            local_comp!(gs, src)
        end
    end

    return gs
end

"""
Perform local complementation on pure graph
"""
function local_comp!(g, id::Int)
    # inneighbors or outneighbors shouldn't really matter because the graph is undirected
    nghbr = inneighbors(g, id)
    for i1 in eachindex(nghbr)
        for i2 in i1+1:length(nghbr)
            toggle_edge!(g, nghbr[i1], nghbr[i2])
        end
    end
    return g
end

"""
Perform local complementation about index id while maintaining the represented quantum state unchanged.
"""
function local_comp!(gs::GraphState{G, SingleQubitOperator} where G, id::Int)
    local_comp!(graph(gs), id)
    # Tweak the VOPs to make the resulting quantum state equivalent
    # See Theorem 3 and Corollary 1 in the paper
    nghbr = inneighbors(graph(gs), id)
    for n in nghbr
        _apply_vop_right!(vops(gs), sPhase(n))
    end
    _apply_vop_right!(vops(gs), sInvSQRTX(id))
    return gs
end

"""
A helper function that applies `ISOLATED_CPHASE_TABLE`.

You should make sure it's suitable to directly copy the result from lookup table. This is the case if either
- vertices are completely isolated
- VOPs of vertices which have non-operand neighbors are all in Z_COMMUTATION_SUBGROUP

An example of the second case is

U₁ = sHadamard, U₂ = sInvPhase

and q2 is connected to some other vertex while q1 has no neighbors except potentially q2.
"""
function apply_isolated_cphase_lookup_table!(gs::GraphState, q1::Int, q2::Int)
    g = graph(gs)

    vs = vops(gs)
    r = ISOLATED_CPHASE_TABLE[(has_edge(g, q1, q2), vs[q1], vs[q2])]
    vs[q1], vs[q2] = vops(r)[1], vops(r)[2]
    # set the edge according to lookup table r
    if has_edge(graph(r), 1, 2) & !has_edge(g, q1, q2)
        add_edge!(g, q1, q2)
    end
    if !has_edge(graph(r), 1, 2) & has_edge(g, q1, q2)
        rem_edge!(g, q1, q2)
    end
end

function apply!(gs::GraphState, gate::sCPHASE)
    # as a convention, we call U₁, U₂ the VOPs of q1, q2 respectively
    q1, q2 = affectedqubits(gate)
    g = graph(gs)

    if has_edge(g, q1, q2)
        # If q1 and q2 are connected, use q2 as a swapping partner of q1 to transfer VOPs of q1 to q2
        remove_vop!(gs, q1, q2)

        # now try remove the VOP of q2 as well

        # find, after the removal of q1's VOP, if q2 has any non-operand neighbors
        nghbr = inneighbors(g, q2)
        targets = filter(x -> x != q1, nghbr)

        if !isempty(targets)
            # If there is at least one non-operand neighbor, remove q2's VOP as well.

            # It's guaranteed that after this step, U₁ ∈ Z_COMMUTATION_SUBGROUP.
            # Because `local_comp!` will only introduce sInvPhase ∈ Z_COMMUTATION_SUBGROUP to non-src and non-target neighbors.
            remove_vop!(gs, q2, targets[1])

            # in which case we are sure both U₁, U₂ ∈ Z_COMMUTATION_SUBGROUP,
            # and we could just toggle the edge between them
            toggle_edge!(g, q1, q2)
        else
            # If we don't have a non-operand swapping partner for q2, the VOP of q1 is identity, and q2 is not connected to other qubits,
            # so we could use the lookup table (see the paper for why even q1's potential non-operand neighbor doesn't matter)
            apply_isolated_cphase_lookup_table!(gs, q1, q2)
        end
   else
        # try reduce both U₁, U₂ to identity
        for q in [q1,q2]
            nghbr = inneighbors(g, q)
            if !isempty(nghbr)
                remove_vop!(gs, q, nghbr[1])
            end
        end

        # now both U₁, U₂ ∈ Z_COMMUTATION_SUBGROUP, we can apply the table
        apply_isolated_cphase_lookup_table!(gs, q1, q2)
   end
   return gs
end
