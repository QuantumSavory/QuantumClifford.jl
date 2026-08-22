@testitem "CSS deterministic ownership" tags=[:qec_determinism] begin
    using Random
    using Test
    using QECCore

    struct CoordinatedCSS <: AbstractCSSCode
        Hx::Matrix{Bool}
        Hz::Matrix{Bool}
        calls::Base.RefValue{Int}
    end

    function QECCore.parity_matrix_xz(c::CoordinatedCSS)
        c.calls[] += 1
        return copy(c.Hx), copy(c.Hz)
    end

    @testset "CSS" begin
        h = QECCore._steane_mat()
        c = CSS(h, h)
        expected_h = copy(h)
        expected_parity = parity_matrix(c)

        @test expected_parity == parity_matrix(Steane7())
        @test parity_matrix_x(c) == expected_h
        @test parity_matrix_z(c) == expected_h
        @test parity_matrix_x(c) isa Matrix{Bool}
        @test parity_matrix_z(c) isa Matrix{Bool}

        h .= false
        @test parity_matrix(c) == expected_parity

        returned_Hx, returned_Hz = parity_matrix_xz(c)
        returned_Hx .= false
        returned_Hz .= false
        @test parity_matrix_x(c) == expected_h
        @test parity_matrix_z(c) == expected_h

        reconstructed = CSS(copy(expected_h), copy(expected_h))
        @test parity_matrix(reconstructed) == expected_parity
        rand(Random.default_rng(), UInt, 100)
        @test parity_matrix(c) == expected_parity

        Hx, Hz = parity_matrix_xz(c)
        @test all(iszero, mod.(Int.(Hx) * transpose(Int.(Hz)), 2))

        coordinated = CoordinatedCSS(copy(expected_h), copy(expected_h), Ref(0))
        @test parity_matrix(coordinated) == expected_parity
        @test coordinated.calls[] == 1
    end
end
