using QuantumClifford
using JET
using ArrayInterface
using Static
using Graphs
using LinearAlgebra

using JET: ReportPass, BasicPass, InferenceErrorReport, UncaughtExceptionReport

# Custom report pass that ignores `UncaughtExceptionReport`
# Too coarse currently, but it serves to ignore the various
# "may throw" messages for runtime errors we raise on purpose
# (mostly on malformed user input)
struct MayThrowIsOk <: ReportPass end

# ignores `UncaughtExceptionReport` analyzed by `JETAnalyzer`
(::MayThrowIsOk)(::Type{UncaughtExceptionReport}, @nospecialize(_...)) = return

# forward to `BasicPass` for everything else
function (::MayThrowIsOk)(report_type::Type{<:InferenceErrorReport}, @nospecialize(args...))
    BasicPass()(report_type, args...)
end

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
        report_pass=MayThrowIsOk(),
        ignored_modules=(
            AnyFrameModule(Graphs.LinAlg),
            AnyFrameModule(Graphs.SimpleGraphs),
            AnyFrameModule(ArrayInterface),
            AnyFrameModule(Static),
            AnyFrameModule(LinearAlgebra)
        )
    )
    @show rep
    @test_broken length(JET.get_reports(rep)) == 0 # false positive from https://github.com/aviatesk/JET.jl/issues/444
    @test length(JET.get_reports(rep)) <= 2
end
