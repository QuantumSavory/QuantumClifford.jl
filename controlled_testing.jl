include("../QuantumClifford/src/QuantumClifford.jl")
using Test
using .QuantumClifford: random_stabilizer, random_destabilizer, @P_str,  @S_str, tCNOT, sCNOT, apply!, sDCNOT, random_pauli, sSWAP, sXNOR, sX, X, sZCX, sZCZ, sHadamard, sInvPhase, tHadamard, sZCY, sPhase, tPhase, sXCY, sYCZ, sYCX, sYCY, sXCX, sXCY, sXCZ

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
    return cat(transform1[control], transform2[target], dims=(1)), cat(inverse2[target] , inverse1[control], dims=(1))
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
    for i in 1:10
        for control in (:X, :Y, :Z)
            for target in (:X, :Y, :Z)
                random_state = random_destabilizer(2)
                implemented_gate = eval(Symbol(:s,control,:C,target))(1,2)

                tr, inv = get_gates(control, target)
                # print(tr, inv)
                test_gates = cat(cat(tr, sCNOT(1,2), dims=(1)), inv, dims=(1))
                @test apply!(copy(random_state), implemented_gate) == walk(test_gates, copy(random_state))
            end
        end
    end
end;