using JET
using ArrayInterface
using Static

function test_jet()
    @testset "JET checks" begin
        @test isempty(JET.get_reports(@report_call random_destabilizer(10)))
        @test isempty(JET.get_reports(@report_call random_stabilizer(10)))
        @test isempty(JET.get_reports(@report_call random_clifford(10)))

        s = random_stabilizer(10)
        p = random_pauli(10)
        c = random_clifford(10)
        c2 = random_clifford(2)
        @test isempty(JET.get_reports(@report_call project!(s,p)))
        @test isempty(JET.get_reports(@report_call apply!(s,p)))
        @test isempty(JET.get_reports(@report_call apply!(s,c)))
        @test isempty(JET.get_reports(@report_call apply!(s,sCNOT(1,2))))
        @test isempty(JET.get_reports(@report_call apply!(s,sHadamard(1))))
        @test isempty(JET.get_reports(@report_call apply!(s,c2,[1,2])))
        @test isempty(JET.get_reports(@report_call canonicalize!(s)))
        @test isempty(JET.get_reports(@report_call canonicalize_gott!(s)))
        @test isempty(JET.get_reports(@report_call canonicalize_rref!(s)))

        @test isempty(JET.get_reports(@report_call QuantumClifford._project!(s,p)))
        @test isempty(JET.get_reports(@report_call QuantumClifford._apply!(s,p)))
        @test isempty(JET.get_reports(@report_call QuantumClifford._apply!(s,c)))
        @test isempty(JET.get_reports(@report_call QuantumClifford._apply!(s,sCNOT(1,2))))
        @test isempty(JET.get_reports(@report_call QuantumClifford._apply!(s,sHadamard(1))))
        @test isempty(JET.get_reports(@report_call QuantumClifford._apply!(s,c2,[1,2])))
        @test isempty(JET.get_reports(@report_call QuantumClifford._canonicalize!(s)))
        @test isempty(JET.get_reports(@report_call QuantumClifford._canonicalize_gott!(s)))
        @test isempty(JET.get_reports(@report_call QuantumClifford._canonicalize_rref!(s,[1,3])))

        rep = report_package("QuantumClifford";
            ignored_modules=(
                AnyFrameModule(Graphs.LinAlg),
                AnyFrameModule(Graphs.SimpleGraphs),
                AnyFrameModule(ArrayInterface),
                AnyFrameModule(Static),
            )
        )
        @show rep
        @test length(JET.get_reports(rep)) == 3
        #= TODO These false positives appear. Figure out how to filter them out.
        ┌ @ /home/stefan/Documents/ScratchSpace/clifford/QuantumClifford/src/linalg.jl:72 LinearAlgebra.rank(::QuantumClifford.Stabilizer)
        │ may throw: QuantumClifford.throw(BadDataStructure("Using a `Stabilizer` type does not ...
        └─────────────────────────────────────────────────────────────────────────────────
        ┌ @ /home/stefan/Documents/ScratchSpace/clifford/QuantumClifford/src/linalg.jl:74 LinearAlgebra.rank(::QuantumClifford.Destabilizer)
        │ may throw: QuantumClifford.throw(BadDataStructure("Using a `Destabilizer` type does not ...
        └─────────────────────────────────────────────────────────────────────────────────
        ┌ @ /home/stefan/Documents/ScratchSpace/clifford/QuantumClifford/src/experimental/NoisyCircuits.jl:387 Base.kwerr(_2, _3, state, noise, indices)
        │┌ @ error.jl:163 Base.kwerr(::Any, ::typeof(QuantumClifford.Experimental.NoisyCircuits.applynoise_branches), ::Register, ::Any, ::Any)
        ││ may throw: Base.throw(Base.MethodError(getproperty(getproperty(getproperty(Base.typeof(applynoise_branches...
        │└────────────────
        =#
    end
end

test_jet()
