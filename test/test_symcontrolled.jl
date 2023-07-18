using Test
using QuantumClifford
using QuantumOptics

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

@testset "Control-Target swap" begin
    for control in (:X, :Z, :Y)
        for target in (:X, :Z, :Y)
            forwgate = eval(Symbol(:s,control,:C,target))(1,2)
            backgate = eval(Symbol(:s,target,:C,control))(1,2)
            forwgatedense = CliffordOperator(forwgate, 2)
            backgatedense = CliffordOperator(backgate, 2)
            @test forwgatedense == tSWAP*backgatedense*tSWAP
        end
    end

    for control in (:X, :Z, :Y)
        for target in (:X, :Z, :Y)
            forwgate = eval(Symbol(:s,control,:C,target))(1,2)
            backgate = eval(Symbol(:s,target,:C,control))(2,1)
            forwgatedense = CliffordOperator(forwgate, 2)
            backgatedense = CliffordOperator(backgate, 2)
            @test forwgatedense == backgatedense
        end
    end

    for (gate1,gate2) in (
        (sCNOT(1,2), sZCX(1,2)),
        (sCPHASE(1,2), sZCZ(1,2)),
    )
        gate1dense = CliffordOperator(gate1, 2)
        gate2dense = CliffordOperator(gate2, 2)
        @test gate1dense == gate2dense
    end
end

@testset "Ket-based definition" begin
    for control in (:X, :Y, :Z)
        for target in (:X, :Y, :Z)
            s = Stabilizer(QuantumClifford._T_str(string(control)))
            k1 = Ket(s)
            s.tab.phases[1] = 0x2
            k2 = Ket(s)
            i = Operator(tId1)
            o = Operator(CliffordOperator(eval(Symbol(:s,target,))(1),1))
            gate = projector(k1)⊗i + (target==:Y ? -im : 1) * projector(k2)⊗o
            implemented_gate = Operator(CliffordOperator(eval(Symbol(:s,control,:C,target))(1,2),2))
            @test gate≈implemented_gate

            target, control = control, target
            s = Stabilizer(QuantumClifford._T_str(string(control)))
            k1 = Ket(s)
            s.tab.phases[1] = 0x2
            k2 = Ket(s)
            i = Operator(tId1)
            o = Operator(CliffordOperator(eval(Symbol(:s,target,))(1),1))
            gate_perm = projector(k1)⊗i + (target==:Y ? -im : 1) * projector(k2)⊗o
            implemented_gate_perm = Operator(CliffordOperator(eval(Symbol(:s,control,:C,target))(1,2),2))
            @test gate_perm≈implemented_gate_perm

            @test permutesystems(gate_perm,[2,1])≈gate
        end
    end
end
