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
end