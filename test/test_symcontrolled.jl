using Test
using QuantumClifford

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

@testset "Controlled gates" begin
    for i in 1:10
        for control in (:X, :Y, :Z)
            for target in (:X, :Y, :Z)
                random_state = random_destabilizer(2)
                implemented_gate = eval(Symbol(:s,control,:C,target))(1,2)

                tr, inv = get_gates(control, target)
                test_gates = [tr..., sCNOT(1,2), inv...]
                @test apply!(copy(random_state), implemented_gate) == mctrajectory!(copy(random_state), test_gates)[1]
            end
        end
    end
end

transforms = Dict(:X => transform_Xbasis(1), :Z => transform_Zbasis(1))

@testset "Change of basis" begin
    for to_basis in (:X, :Y, :Z)
        for from_basis in (:X, :Z)
            start_state = Stabilizer(QuantumClifford._T_str(string(from_basis)))
            forward, backward = transforms[from_basis][1][to_basis], transforms[from_basis][2][to_basis]
            mid_state = mctrajectory!(copy(start_state), backward)[1]
            end_state = mctrajectory!(copy(mid_state), forward)[1]
            @test start_state == end_state
            @test mid_state == Stabilizer(QuantumClifford._T_str(string(to_basis)))
        end
    end
end

@testset "Explicit control" begin
    for control in (:X, :Z, :Y)
        for target in (:X, :Z, :Y)
            for targetstate in ("X","Y","Z")
                twoqgate = eval(Symbol(:s,control,:C,target))(1,2)
                oneqgate = eval(Symbol(:s,target))(2)
                state0str = string(control)*"I I"*targetstate
                state1str = "-"*string(control)*"I I"*targetstate
                state0 = Stabilizer(QuantumClifford._T_str(string(state0str)))
                state1 = Stabilizer(QuantumClifford._T_str(string(state1str)))
                @test canonicalize!(twoqgate*state0) == canonicalize!(state0)
                @test canonicalize!(twoqgate*state1) == canonicalize!(oneqgate*state1)
            end
        end
    end
end