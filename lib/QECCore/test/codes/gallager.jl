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
                        c = GallagerLDPC(μ, wc, wr)
                        H = parity_matrix(c)
                        c1 = GallagerLDPC(μ, wc, wr)
                        H1 = parity_matrix(c1)
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
                        c = GallagerLDPC(μ, wc, wr, 1)
                        H = parity_matrix(c)
                        c1 = GallagerLDPC(μ, wc, wr, 2)
                        H1 = parity_matrix(c1)
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
