using JET
using Base.Threads: nthreads

function test_jet()
    @testset "JET checks" begin
        n = nthreads()

        @test_broken isempty(JET.get_reports(@report_call random_destabilizer(10)))
        @test_broken isempty(JET.get_reports(@report_call random_stabilizer(10)))
        @test_broken isempty(JET.get_reports(@report_call random_clifford(10)))

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

        @test_broken isempty(JET.get_reports(report_package("QuantumClifford"; ignored_modules=(AnyFrameModule(Graphs.LinAlg),))))
    end
end

test_jet()
