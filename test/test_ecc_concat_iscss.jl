@testitem "Concatenated CSS parity matrices" tags=[:ecc, :ecc_base] begin
    using QuantumClifford.ECC: Concat, CSS, Shor9, Steane7, Toric,
        code_n, iscss, parity_matrix, parity_matrix_x, parity_matrix_z

    @test iscss(Concat(Steane7(), Steane7())) === true

    cases = (
        (Concat(Shor9(), Steane7()), (17, 63), (45, 63)),
        (Concat(Toric(2, 2), Shor9()), (31, 72), (39, 72)),
        (Concat(CSS(falses(0, 2), Bool[1 1]), Steane7()), (3, 14), (10, 14)),
    )
    for (code, xsize, zsize) in cases
        Hx = parity_matrix_x(code)
        Hz = parity_matrix_z(code)
        H = parity_matrix(code)
        n = code_n(code)
        xpart = @view H[:, 1:n]
        zpart = @view H[:, n+1:end]
        xrows = vec(any(xpart; dims=2))
        zrows = vec(any(zpart; dims=2))

        @test size(Hx) == xsize
        @test size(Hz) == zsize
        @test eltype(Hx) === Bool
        @test eltype(Hz) === Bool
        @test all(xrows .⊻ zrows)
        @test Hx == xpart[xrows, :]
        @test Hz == zpart[zrows, :]
    end
end
