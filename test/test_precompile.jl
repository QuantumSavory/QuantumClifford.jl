@testitem "Precompile" begin
    using Random

    rng = Random.default_rng()
    saved_rng = copy(rng)
    try
        QuantumClifford._precompile_()
        expected_rng = copy(saved_rng)
        @test rand(rng, UInt64) == rand(expected_rng, UInt64)
    finally
        copy!(rng, saved_rng)
    end
end
