using Random
using QuantumClifford

using QuantumClifford: stab_looks_good, destab_looks_good, mixed_stab_looks_good, mixed_destab_looks_good
using QuantumClifford: apply_single_x!, apply_single_y!, apply_single_z!
using InteractiveUtils
using Test

test_size = 10

function transform_Zbasis(qubit)
    transformations = Dict(:X => [sHadamard(qubit),], :Y => [sInvPhase(qubit),sHadamard(qubit)], :Z => [sHadamard(qubit), sHadamard(qubit)])
    inverse_transformations = Dict(:X => [sHadamard(qubit),], :Y => [sHadamard(qubit), sPhase(qubit)], :Z => [sHadamard(qubit), sHadamard(qubit)])
    return transformations, inverse_transformations
end

function transform_Xbasis(qubit)
    transformations = Dict(:Z => [sHadamard(qubit),], :Y => [sInvPhase(qubit),], :X => [sHadamard(qubit), sHadamard(qubit)])
    inverse_transformations = Dict(:Z => [sHadamard(qubit)], :Y => [sPhase(qubit),], :X => [sHadamard(qubit), sHadamard(qubit)])
    return transformations, inverse_transformations
end

function get_gates(control, target)
    transform1, inverse1 = transform_Zbasis(1) # get gates which swith control from Z basis to specified
    transform2, inverse2 = transform_Xbasis(2) # get gates which swith control from X basis to specified
    return [transform1[control]..., transform2[target]...], [inverse2[target]..., inverse1[control]...]
end

function walk(gates, state)
    cstate = copy(state)
    for gate in gates
        cstate = apply!(cstate, gate)
    end
    return cstate
end

@testset "Controlled gates" begin
    # test ten time for each
    for i in 1:test_size
        for control in (:X, :Y, :Z)
            for target in (:X, :Y, :Z)
                random_state = random_destabilizer(2)
                implemented_gate = eval(Symbol(:s,control,:C,target))(1,2)

                tr, inv = get_gates(control, target)
                # print(tr, inv)
                test_gates = [tr..., sCNOT(1,2), inv...]
                @test apply!(copy(random_state), implemented_gate) == walk(test_gates, copy(random_state))
            end
        end
    end
end

transforms = Dict(:X => transform_Xbasis(1), :Z => transform_Zbasis(1))

@testset "Change of basis" begin
    for from_basis in (:X, :Z)
        for to_basis in (:X, :Y, :Z)
            start_state = Stabilizer(QuantumClifford._T_str(string(from_basis)))
            forward, backward = [t[to_basis] for t in transforms[from_basis]]
            mid_state = walk(forward, start_state)
            end_state = walk(backward, mid_state)
            @test start_state == end_state
            @test mid_state == Stabilizer(QuantumClifford._T_str(string(to_basis)))
        end
    end
end