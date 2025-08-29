@testitem "ECC Trivariate Tricycle Code" tags=[:ecc] begin
    @static if !Sys.iswindows() && Sys.ARCH == :x86_64 && VERSION >= v"1.11"
        using Oscar
        using QECCore
        import HiGHS
        import JuMP
        using Nemo: matrix, GF
        using QuantumClifford: stab_looks_good, stab_to_gf2
        using QuantumClifford.ECC

        @testset "Trivariate Tricycle Codes" begin
            F2 = GF(2)
            R, (x, y, z) = polynomial_ring(F2, [:x, :y, :z])
            table_i = [
                # n   , k  ,  ℓ, m, p, A                           , B                            ,   C
                (72   , 6  ,  4, 3, 2, 1 + y         + x*y^2       , 1 + y*z       + x^2*y^2      , 1 + x*y^2*z     + x^2*y      ),
                (180  , 12 ,  5, 4, 3, 1 + x^2*y^3*z + x^4*y       , 1 + x^3       + x^4*z^2      , 1 + x^3*y^3     + x^4*y*z^2  ),
                (432  , 12 ,  6, 6, 4, 1 + x*y*z^3   + x^3*y^4*z^2 , 1 + x^3*y*z^2 + x^3*y^2*z^3  , 1 + x^4*y^3*z^3 + x^5*z^2    ),
            ]

            table_iii = [
                # n , k , ℓ, m, p, A         , B        ,   C
                (36 , 3 , 3, 2, 2, 1 + x*y*z , 1 + x^2*z, 1 + x^2*y    ),
                (48 , 3 , 4, 2, 2, 1 + x     , 1 + x*z  , 1 + x*y      ),
                (54 , 3 , 3, 3, 2, 1 + y*z   , 1 + x*z  , 1 + x*y*z    ),
                (60 , 3 , 5, 2, 2, 1 + x*z   , 1 + x*y  , 1 + x*y*z    ),
                (90 , 3 , 5, 3, 2, 1 + x     , 1 + x*y  , 1 + x^2*y^2*z),
            ]

            table_iv = [
                # n  , k  , ℓ, m, p, A                             , B             , C
                (36  , 6  , 3, 2, 2, (1 +       z)*(1 + x        ) , 1 + x         , 1 + x*y*z       ),
                (48  , 6  , 4, 2, 2, (1 + x^2*y*z)*(1 + x*z      ) , 1 + x^3       , 1 + x^3*y*z     ),
                (54  , 9  , 3, 3, 2, (1 +       z)*(1 + y        ) , 1 + x*y*z     , 1 + x^2*y^2     ),
                (108 , 6  , 4, 3, 3, (1 + x^2    )*(1 + x*z      ) , 1 + x^2*y^2   , 1 + x^2*y^2*z^2 ),
                (108 , 18 , 4, 3, 3, (1 + x^2    )*(1 + x*y^2*z^2) , 1 + x^2*y^2*z , 1 + x^2*y^2*z   ),
            ]

            table_v = [
                # n  , k ,  ℓ, m, p, A                          , B                          ,   C
                (36  , 6 ,  3, 2, 2, 1 + x         + x^2*z      , 1 + x*y       + x^2*y      , 1 + x*y*z       + x^2        ),
                (72  , 6 ,  4, 3, 2, 1 + y         + x*y^2      , 1 + y*z       + x^2*y^2    , 1 + x*y^2*z     + x^2*y      ),
                (81  , 6 ,  3, 3, 3, 1 + x         + x*y        , 1 + y         + y*z        , 1 + z           + z*x        ),
                (126 , 6 ,  7, 3, 2, 1 + y^2       + x^4*y      , 1 + x*y       + x^4*y^2*z  , 1 + x^2*y*z     + x^2*y^2    ),
                (135 , 12,  5, 3, 3, 1 + x         + x^3*y^2    , 1 + x*z^2     + x^2*y*z    , 1 + x*y^2*z     + x^2*y*z    ),
                (180 , 12,  5, 4, 3, 1 + x^2*y^3*z + x^4*y      , 1 + x^3       + x^4*z^2    , 1 + x^3*y^3     + x^4*y*z^2  ),
                (288 , 6 ,  8, 4, 3, 1 + y*z       + x^7*y*z^2  , 1 + x^2*y*z   + x^5*z^2    , 1 + x^2*y^3*z   + x^6*y^2*z^2),
                (432 , 12,  6, 6, 4, 1 + x*y*z^3   + x^3*y^4*z^2, 1 + x^3*y*z^2 + x^3*y^2*z^3, 1 + x^4*y^3*z^3 + x^5*z^2    ),
                (588 , 9 ,  7, 7, 4, 1 + y*z       + x^3*y^2*z  , 1 + x*y^4*z   + x^4*y^4*z^2, 1 + x^2*y^2*z   + x^2*y^4    ),
                (648 , 6 ,  6, 6, 6, 1 + x*z^4     + x^3*y*z^4  , 1 + x*y^4*z   + x^5*y*z^5  , 1 + x^2*y*z^2   + x^2*y^2*z^3),
                (648 , 12,  6, 6, 6, 1 + x         + x^4*z      , 1 + y         + x*y*z^4    , 1 + z           + x^3*y^2*z  ),
                (735 , 18,  7, 7, 5, 1 + y^2       + x^6*y^5    , 1 + x^2*y*z   + x^2*y^4*z^3, 1 + x^4*y^5*z^2 + x^5*y      ),
                (840 , 9 ,  8, 7, 5, 1 + y^3*z     + x*y*z^2    , 1 + y^5*z     + x^3*y^4*z  , 1 + x*y^3*z^4   + x^7*y*z    ),
                (1029, 18,  7, 7, 7, 1 + y^2*z^3   + x^5*y^5*z  , 1 + x^2*y^5*z + x^4*z^5    , 1 + x^4*y^6     + x^5*y^4*z^6),
            ]

                table_vi = [
                # n   , k  ,  ℓ, m, p, A                       , B                         ,   C
                (108  , 18 ,  4, 3, 3, (1 + x^2)*(1 + x*y*z^2) , (1 + x^2)*(1 + x^3*y^2*z) ,  1  + x^2*y^2*z^2      ),
                (108  , 60 ,  4, 3, 3, (1 + x^2)*(1 +       z) , (1 + x^2)*(1 + x^2*y*z^2) , (1  + x^2)*(1 + z^2)   )
            ]

            for (n, k, ℓ, m, p, A_poly, B_poly, C_poly) in vcat(table_i, table_iii, table_iv, table_v, table_vi)
                F2 = GF(2)
                R, (x, y, z) = polynomial_ring(F2, [:x, :y, :z])
                I = ideal(R, [x^ℓ - 1, y^m - 1, z^p - 1])
                S, _ = quo(R, I)
                A = S(A_poly)
                B = S(B_poly)
                C = S(C_poly)
                c = TrivariateTricycleCode(ℓ, m, p, A, B, C)
                stab = parity_checks(c)
                mat = matrix(GF(2), stab_to_gf2(stab))
                computed_rank = rank(mat)
                @test computed_rank == code_n(c) - code_k(c)
                # A TT code is defined on n = 3*ℓ*m*p data qubits.
                @test code_n(c) == n == code_n(stab) == 3*ℓ*m*p
                @test code_k(c) == k == code_k(stab)
                @test stab_looks_good(stab, remove_redundant_rows=true) == true
                @test iszero(mod.(metacheck_matrix_z(c)*parity_matrix_z(c), 2))
            end
        end
    end
end
