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
        if n>1
            @test_broken isempty(JET.get_reports(@report_call apply!(s,c)))
        else
            @test isempty(JET.get_reports(@report_call apply!(s,c)))
        end
        @test_broken isempty(JET.get_reports(@report_opt apply!(s,sCNOT(1,2))))
        @test_broken isempty(JET.get_reports(@report_opt apply!(s,sHadamard(1))))
        @test_broken isempty(JET.get_reports(@report_opt apply!(s,c2,[1,2])))
        @test isempty(JET.get_reports(@report_call canonicalize!(s)))
        @test isempty(JET.get_reports(@report_call canonicalize_gott!(s)))
        @test isempty(JET.get_reports(@report_call canonicalize_rref!(s)))

        @test isempty(JET.get_reports(@report_opt QuantumClifford._project!(s,p)))
        @test isempty(JET.get_reports(@report_opt QuantumClifford._apply!(s,p)))
        @test_broken isempty(JET.get_reports(@report_opt QuantumClifford._apply!(s,c)))
        if n>1
            @test_broken isempty(JET.get_reports(@report_opt QuantumClifford._apply!(s,sCNOT(1,2))))
            @test_broken isempty(JET.get_reports(@report_opt QuantumClifford._apply!(s,sHadamard(1))))
        else
            @test isempty(JET.get_reports(@report_opt QuantumClifford._apply!(s,sCNOT(1,2))))
            @test isempty(JET.get_reports(@report_opt QuantumClifford._apply!(s,sHadamard(1))))
        end
        @test_broken isempty(JET.get_reports(@report_opt QuantumClifford._apply!(s,c2,[1,2])))
        @test isempty(JET.get_reports(@report_opt QuantumClifford._canonicalize!(s)))
        @test isempty(JET.get_reports(@report_opt QuantumClifford._canonicalize_gott!(s)))
        @test isempty(JET.get_reports(@report_opt QuantumClifford._canonicalize_rref!(s,[1,3])))

        @test_broken isempty(JET.get_reports(report_package("QuantumClifford")))
    end
end

test_jet()
