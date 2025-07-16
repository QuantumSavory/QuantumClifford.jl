@testitem "Codes" begin
    using QECCore
    using Test
    using TestItemRunner
    
    @testset "fivequbit" begin
        include("codes/fivequbit.jl")
    end

    @testset "repetition" begin
        include("codes/repetition.jl")
    end

    @testset "toric" begin
        include("codes/toric.jl")
    end

    @testset "surface" begin
        include("codes/surfacecode.jl")
    end

    @testset "steanecode" begin
        include("codes/steanecode.jl")
    end

    @testset "shorcode" begin
        include("codes/shorcode.jl")
    end

    @testset "clevecode" begin
        include("codes/clevecode.jl")
    end

    @testset "css" begin
        include("codes/css.jl")
    end

    @testset "quantumreedmuller" begin
        include("codes/quantumreedmuller.jl")
    end

    @testset "reedmuller" begin
        include("codes/reedmuller.jl")
    end

    @testset "color_codes" begin
        include("codes/color_codes.jl")
    end

    @testset "gottesman" begin
        include("codes/gottesman.jl")
    end

    @testset "golay" begin
        include("codes/golay.jl")
    end

    @testset "hamming" begin
        include("codes/hamming.jl")
    end

    @testset "Delfosse-Reichardt Repetition Code" begin
        include("codes/delfosse_reichardt_repcode.jl")
    end

    @testset "Delfosse-Reichardt Code" begin
        include("codes/delfosse_reichardt_code.jl")
    end

    @testset "Delfosse-Reichardt Generalized [[8,2,3]] Code" begin
        include("codes/delfosse_reichardt_823_code.jl")
    end
end
