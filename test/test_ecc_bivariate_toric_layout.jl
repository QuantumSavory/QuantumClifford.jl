@testitem "ECC Bivaraite Toric Layout" begin
    using QuantumClifford.ECC
    using QuantumClifford.ECC: AbstractECC, bivariate_toric_layout

    codes = [[6 ,6 ,3 ,1 ,2 ,3 ,1 ,2 ],
             [15,3 ,9 ,1 ,2 ,0 ,2 ,7 ],
             [6 ,9 ,3 ,1 ,2 ,3 ,1 ,2 ],
             [12,6 ,3 ,1 ,2 ,3 ,1 ,2 ],
             [12,12,3 ,2 ,7 ,3 ,1 ,2 ],
             [30,6 ,9 ,1 ,2 ,3 ,25,26],
             [21,18,3 ,10,17,5 ,3 ,19]]
    # cross checks taken from https://arxiv.org/pdf/2404.17676.
    @test bivariate_toric_layout(codes[1]) == [(6, 6, 0, 2), (6, 6, 0, 3), (6, 6, 0, 4), (6, 6, 0, 5), (6, 6, 1, 2), (6, 6, 1, 3), (6, 6, 1, 4), (6, 6, 1, 5), (6, 6, 2, 0), (6, 6, 2, 1), (6, 6, 2, 2), (6, 6, 2, 3), (6, 6, 3, 0), (6, 6, 3, 1), (6, 6, 3, 2), (6, 6, 3, 3), (6, 6, 4, 0), (6, 6, 4, 1), (6, 6, 4, 4), (6, 6, 4, 5), (6, 6, 5, 0), (6, 6, 5, 1), (6, 6, 5, 4), (6, 6, 5, 5)]
    @test  bivariate_toric_layout(codes[2]) == [(15, 3, 0, 2), (15, 3, 0, 3), (15, 3, 1, 2), (15, 3, 1, 3), (3, 15, 2, 0), (3, 15, 2, 1), (3, 15, 2, 4), (3, 15, 2, 5), (3, 15, 3, 0), (3, 15, 3, 1), (3, 15, 3, 4), (3, 15, 3, 5), (15, 3, 4, 2), (15, 3, 4, 3), (15, 3, 5, 2), (15, 3, 5, 3)]
    @test bivariate_toric_layout(codes[3]) == [(18, 3, 0, 4), (18, 3, 0, 5), (18, 3, 1, 4), (18, 3, 1, 5), (9, 6, 2, 0), (9, 6, 2, 1), (9, 6, 2, 2), (9, 6, 2, 3), (9, 6, 3, 0), (9, 6, 3, 1), (9, 6, 3, 2), (9, 6, 3, 3), (18, 3, 4, 4), (18, 3, 4, 5), (18, 3, 5, 4), (18, 3, 5, 5)]
    @test bivariate_toric_layout(codes[4]) == [(12, 6, 0, 4), (12, 6, 0, 5), (12, 6, 1, 4), (12, 6, 1, 5), (6, 12, 2, 0), (6, 12, 2, 1), (6, 12, 2, 2), (6, 12, 2, 3), (6, 12, 3, 0), (6, 12, 3, 1), (6, 12, 3, 2), (6, 12, 3, 3), (12, 6, 4, 4), (12, 6, 4, 5), (12, 6, 5, 4), (12, 6, 5, 5)]
    @test bivariate_toric_layout(codes[5]) == [(12, 12, 0, 0), (12, 12, 0, 1), (12, 12, 0, 4), (12, 12, 0, 5), (12, 12, 1, 0), (12, 12, 1, 1), (12, 12, 1, 4), (12, 12, 1, 5), (12, 12, 2, 0), (12, 12, 2, 1), (12, 12, 2, 2), (12, 12, 2, 3), (12, 12, 3, 0), (12, 12, 3, 1), (12, 12, 3, 2), (12, 12, 3, 3), (12, 12, 4, 2), (12, 12, 4, 3), (12, 12, 4, 4), (12, 12, 4, 5), (12, 12, 5, 2), (12, 12, 5, 3), (12, 12, 5, 4), (12, 12, 5, 5)]
    @test bivariate_toric_layout(codes[6]) == [(6, 30, 2, 2), (6, 30, 2, 3), (6, 30, 3, 2), (6, 30, 3, 3), (30, 6, 4, 0), (30, 6, 4, 1), (30, 6, 5, 0), (30, 6, 5, 1)]
    @test bivariate_toric_layout(codes[7]) == [(18, 21, 2, 2), (18, 21, 2, 3), (18, 21, 3, 2), (18, 21, 3, 3)]
end
