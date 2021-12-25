function test_gf2()
    @testset "GF(2) representations" begin
        @testset "Equivalence of GF(2) Gaussian elimination and Stabilizer canonicalization" begin
            for n in test_sizes
                for rep in 1:5
                    s = random_stabilizer(n)[randperm(n)[1:rand(n÷2+1:n)]]
                    cs = canonicalize!(copy(s));
                    H = stab_to_gf2(cs);
                    cH = gf2_gausselim!(stab_to_gf2(s));
                    @test H==cH
                end
            end
        end
        @testset "GF(2) H and G matrices" begin
            for n in test_sizes
                for rep in 1:5
                    H = random_invertible_gf2(n)[randperm(n)[1:rand(n÷2+1:n)],:]
                    H = gf2_gausselim!(H)
                    G = gf2_H_to_G(H)
                    @test sum(G*H' .%2)==0;
                end
            end
        end
    end    
end

test_gf2()