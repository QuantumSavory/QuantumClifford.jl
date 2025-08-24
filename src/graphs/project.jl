# Some documentation for myself:
# In line with the existing definition
# - `projectrand!(gs::GraphState, m::sMX)` doesn't care about the register/classical bit
# - `apply!(state::Register, m::sMX)` does, and it depends on `projectrand!`
# And these should be the only public facing APIs for now

function remove_neighbors!(g, qubit::Int)
    for nghbr in copy(inneighbors(g, qubit))
        rem_edge!(g, qubit, nghbr)
    end
end

"""
Project a pure graph state (all VOPs are identity) onto eigenspaces of the Z observable.
"""
function _projectZ_pure_graph!(state::GraphState, qubit::Int, res::Bool)
    # For a = qubit to measure, b ∈ nghbrs(a)
    # Cₐ -> Cₐ * Xₐ^res * Hₐ
    # Cᵦ -> Cᵦ * Zᵦ^res
    # And remove all edeges between a and its neighbors
    g = graph(state)
    v = vops(state)

    # Tweak VOPs
    if res
        _apply_vop_right!(v, sX(qubit))
    end
    _apply_vop_right!(v, sHadamard(qubit))
    for nghbr in inneighbors(g, qubit)
        if res
            _apply_vop_right!(v, sZ(nghbr))
        end
    end

    # Modify graph
    remove_neighbors!(g, qubit)

    return state
end

"""
Project a pure graph state onto eigenspaces of the Y observable.
"""
function _projectY_pure_graph!(state::GraphState, qubit::Int, res::Bool)
    g = graph(state)
    v = vops(state)

    function tweak_vop(q)
        if res
            _apply_vop_right!(v, sInvPhase(q))
        else
            _apply_vop_right!(v, sPhase(q))
        end
    end

    # Tweak VOPs of target qubit and its neighbors
    tweak_vop(qubit)
    for b in inneighbors(g, qubit)
        tweak_vop(b)
    end

    # perform local complementation on edge
    local_comp!(g, qubit)

    # remove the neighbors of the target qubit, this step is missing in the paper
    remove_neighbors!(g, qubit)

    return state
end

"""
Project a pure graph state onto eigenspaces of X observable.
"""
function _projectX_pure_graph!(state::GraphState, qubit::Int, res::Bool)
    g = graph(state)
    v = vops(state)
    a = qubit

    # Cₐ -> Cₐ Z^res
    if res
        _apply_vop_right!(v, sZ(a))
    end

    nghbrs_a = copy(inneighbors(g, a))
    if !isempty(nghbrs_a)
        b = nghbrs_a[1]
        nghbrs_b = copy(inneighbors(g, b))
        nghbr_a_and_b = intersect(nghbrs_a, nghbrs_b)

        # tweak the VOPs first
        # C_b -> C_b (sInvSQRTY)(†) where † is applied only if res is true
        if res
            _apply_vop_right!(v, sSQRTY(b))
        else
            _apply_vop_right!(v, sInvSQRTY(b))
        end

        if res
            # nghbr(b) \ nghbr(a) \ {a}
            for c in setdiff(nghbrs_b, nghbrs_a, a)
                _apply_vop_right!(v, sZ(c))
            end
        else
            # nghbr(a) \ nghbr(b) \ {b}
            for c in setdiff(nghbrs_a, nghbrs_b, b)
                _apply_vop_right!(v, sZ(c))
            end
        end

        # complementation on edges {{c, d} | c ∈ nghbr(a), d ∈ nghbr(b)}
        for (c, d) in Iterators.product(nghbrs_a, nghbrs_b)
            if (c ∈ nghbrs_b) && (d ∈ nghbrs_a)
                if c < d
                    toggle_edge!(g, c, d)
                end
            else
                toggle_edge!(g, c, d)
            end
        end

        # complementation on edges {{c, d} | c, d ∈ nghbr(a) ∩ nghbr(b)}
        for c in eachindex(nghbr_a_and_b)
            for d in c+1:length(nghbr_a_and_b)
                toggle_edge!(g, nghbr_a_and_b[c], nghbr_a_and_b[d])
            end
        end

        # complement on edges {{b, d} | d ∈ nghbr(a) \ {b}}
        for d in nghbrs_a[2:end]
            toggle_edge!(g, b, d)
        end

        remove_neighbors!(g, a)
    end
    return state
end

function _project_graph!(state::GraphState, qubit::Int, op::Stabilizer, res::Bool)
    # Cₐ† opₐ Cₐ
    op_prime = apply!(op, inv(vops(state)[qubit]))

    ispositive = (phases(op_prime)[1]) != 0x02
    # normalize the phase to match later
    phases(op_prime)[1] = 0x00
    res_tilde = ispositive ? res : !res

    if op_prime == S"Z"
        _projectZ_pure_graph!(state, qubit, res_tilde)
    elseif op_prime == S"X"
        _projectX_pure_graph!(state, qubit, res_tilde)
    elseif op_prime == S"Y"
        _projectY_pure_graph!(state, qubit, res_tilde)
    end

    return state
end