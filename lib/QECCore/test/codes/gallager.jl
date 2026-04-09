@testitem "Gallager LDPC code" begin
    using LinearAlgebra
    using QECCore
    using QECCore: generator
    using Random
    using SparseArrays
    using Nemo: matrix, GF, echelon_form

    @testset "Gallager LDPC" begin
        for μ in 2:10
            for wc in 3:10
                for wr in [wc+1, wc+3, 2*wc]
                    wr <= wc && continue
                    @testset "μ=$μ, wc=$wc, wr=$wr" begin
                        rng1 = MersenneTwister(123)
                        rng2 = MersenneTwister(123)
                        H = random_Gallager_ldpc(rng1, μ, wc, wr)
                        H1 = random_Gallager_ldpc(rng2, μ, wc, wr)
                        @test Matrix(H) == Matrix(H1)
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
                        r = rank(matrix(GF(2), H))
                        @test r ≤ m && 0 < 1 - r/n < 1
                    end
                    @testset "μ=$μ, wc=$wc, wr=$wr" begin
                        rng1 = MersenneTwister(1)
                        rng2 = MersenneTwister(2)
                        H = random_Gallager_ldpc(rng1, μ, wc, wr)
                        H1 = random_Gallager_ldpc(rng2, μ, wc, wr)
                        @test Matrix(H) != Matrix(H1)
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
                        r = rank(matrix(GF(2), H))
                        @test r ≤ m && 0 < 1 - r/n < 1
                    end
                end
            end
        end
    end
end
