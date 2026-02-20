@testitem "ECC ZSZ" tags=[:ecc, :ecc_bespoke_checks, :oscar_required] begin
    using Oscar
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: DistanceMIPAlgorithm, code_n, code_k, distance, ZSZCode
    using Random

    @testset "Simple ZSZ Code Construction" begin
        # Setting a seed so this doesn't flake out on us in CI.
        Random.seed!(42)
        # I'm using the example Z_5 x| Z_4 with q=2.
        # Checked that 2^4 = 16 which is 1 mod 5, so we're good.
        # This gives us a group of size 20.
        l, m, q = 5, 4, 2
        
        found_good_code = false
        
        # Let's give it 20 shots to find something decent.
        for _ in 1:20
            # Helper to just grab unique random elements represented by (val_x, val_y)
            function rand_poly_tuples()
                elts = [(0, 0)] # identity element x^0 * y^0
                while length(elts) < 3
                    # Picking random powers for x and y
                    i = rand(0:l-1)
                    j = rand(0:m-1)
                    if !((i, j) in elts)
                        push!(elts, (i, j))
                    end
                end
                return elts
            end

            # Generating my random polynomials A and B
            A_tups = rand_poly_tuples()
            B_tups = rand_poly_tuples()

            # Constructing the 2BGA code with them
            c = ZSZCode(l, m, q, A_tups, B_tups)
            
            # Now let's check if it actually worked.
            # n should be 2 * size of group (2 * 20 = 40).
            n = code_n(c)
            k = code_k(c)
            
            # I want a code that actually encodes something (k > 0).
            if n == 2*l*m && k > 0
                @test n == 40
                @test k >= 2 
                found_good_code = true
                break
            end
        end
        
        @test found_good_code
    end
    
    @testset "Invalid ZSZ Parameters" begin
        # Just making sure I catch bad parameters.
        # 2^3 = 8 which is 3 mod 5, not 1. Should fail.
        @test_throws ArgumentError ZSZCode(5, 3, 2, [(0,0)], [(0,0)])
    end
end
