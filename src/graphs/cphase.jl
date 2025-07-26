using Graphs: SimpleGraph, inneighbors, has_edge, rem_edge!, add_edge!

function toggle_edge!(g::SimpleGraph, v1::Int, v2::Int)
    if has_edge(g, v1, v2)
        rem_edge!(g, v1, v2)
    else
        add_edge!(g, v1, v2)
    end
    return g
end

# Perform Local Complementation repeatedly to remove the VOP of src to the target
function remove_vop!(gs::GraphState, src::Int, target::Int)
    # Decomposition of the source VOP in the order of right to left
    src_seq = IP_SQRTX_DECOMPOSITION_TABLE[vops(gs)[src]]

    if !has_edge(graph(gs), src, target)
        throw("Target qubit must be a neighbor of the source qubit")
    end

    for v in src_seq
        if v isa sInvPhase
            # perform local_comp on the target qubit
            local_comp!(gs, target)
        elseif v isa sSQRTX
            local_comp!(gs, src)
        end
    end

    return gs
end

# Perform local complementation on graph
function local_comp!(g::SimpleGraph, id::Int)
    # inneighbors or outneighbors shouldn't really matter because the graph is undirected
    nghbr = inneighbors(g, id)
    for i1 in eachindex(nghbr)
        for i2 in i1+1:length(nghbr)
            toggle_edge!(g, nghbr[i1], nghbr[i2])
        end
    end
    return g
end

# Perform Local Complementation about index id, producing an equivalent quantum state
function local_comp!(gs::GraphState{<:SimpleGraph, SingleQubitOperator}, id::Int)
    local_comp!(graph(gs), id)
    # Tweak the VOPs to make the resulting quantum state equivalent
    nghbr = inneighbors(graph(gs), id)
    for n in nghbr
        _apply_vop_right!(vops(gs), sPhase(n))
    end
    _apply_vop_right!(vops(gs), sInvSQRTX(id))
    return gs
end

function apply!(gs::GraphState, gate::sCPHASE)
   q1, q2 = affectedqubits(gate)
   g = graph(gs)

   if has_edge(g, q1, q2)
        # q1 and q2 are connected
        # remove the VOP of q1 to q2
        remove_vop!(gs, q1, q2)
        # try remove the VOP of q2 as well
        nghbr = inneighbors(g, q2)
        targets = filter(x -> x != q1, nghbr)
        if !isempty(targets)
            remove_vop!(gs, q2, targets[1])
            # in which case we are sure both VOPs of q1 and q2 are in the Z_COMMUTATION_SUBGROUP,
            # and we could just toggle the edge between them
            toggle_edge!(g, q1, q2)
        else
            # the VOP of q1 is identity, and q2 is not connected to other qubits,
            # so we could use the lookup table (see the paper for why)
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
   else
        # try remove the VOP of q1 or q2, whoever has the neighbor first
        nghbr1 = inneighbors(g, q1)
        nghbr2 = inneighbors(g, q2)
        nghbr = isempty(nghbr2) ? nghbr1 : nghbr2
        src = isempty(nghbr2) ? q1 : q2
        if !isempty(nghbr)
            # we have soemthing to use to remove
            remove_vop!(gs, src, nghbr[1])
        end

        # in either case we can apply the table
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
   return gs
end
