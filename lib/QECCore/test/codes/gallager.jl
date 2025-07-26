@testitem "Gallager LDPC code" begin
    using LinearAlgebra
    using QECCore
    using QECCore: generator
    using Random
    using SparseArrays
    using Nemo: matrix, GF, echelon_form

    rng = MersenneTwister()
    @testset "Gallager LDPC" begin
        for μ in 1:20
            for wc in 3:20
                for wr in [wc+1, wc+3, 2*wc]
                    wr <= wc && continue
                    @testset "μ=$μ, wc=$wc, wr=$wr" begin
                        H = Matrix(random_gallager_ldpc_code(μ, wc, wr; rng))
                        m, n = size(H)
                        @test m == μ*wc
                        @test n == μ*wr
                        @test all(sum(H, dims=1) .== wc)
                        @test all(sum(H, dims=2) .== wr)
                        H1 = H[1:μ, :]
                        for d in 2:wc
                            block = H[(d-1)*μ+1:d*μ, :]
                            for row in 1:μ
                                @test sort(block[row, :]) == sort(H1[row, :])
                            end
                        end
                        @test all(x -> x == 0 || x == 1, H)
                        mat = matrix(GF(2), H)
                        r = rank(mat)
                        @test r ≤ m
                        rate = 1 - r/n
                        @test 0 < rate < 1
                    end
                end
            end
        end
    end
end
