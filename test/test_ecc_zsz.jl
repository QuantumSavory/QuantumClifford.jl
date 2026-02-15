@testitem "ECC ZSZ" tags=[:ecc, :ecc_bespoke_checks] begin
    using Oscar
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: DistanceMIPAlgorithm, code_n, code_k, distance, twobga_from_fp_group
    using Random

    # So, I'm implementing the ZSZ codes from the Guo, Hong, Lucas (2025) paper.
    # Basically, Section 3 defines this ZSZ group Z_l x| Z_m.
    # The presentation is: < x, y | x^l = 1, y^m = 1, y x y^-1 = x^q > 
    # and we gotta make sure q^m = 1 mod l for it to make sense.

    # Here's a helper I wrote to spin up the Z_l x| Z_m group rings.
    function zsz_group_rings(l, m, q)
        # First things first, checking Eq 8 from the paper: q^m must be 1 mod l.
        # If not, I'm throwing an error cause that's invalid.
        if powermod(q, m, l) != 1
            throw(ArgumentError("Condition q^m = 1 (mod l) not satisfied for l=$l, m=$m, q=$q"))
        end

        # Now I'm building the free group with 2 generators.
        G = free_group(2)
        x, y = gens(G)
        
        # Defining the relations like the paper says.
        # x^l is identity, y^m is identity, and the conjugation rule.
        rels = [x^l, y^m, y*x*y^-1 * x^-q]
        Q, _ = quo(G, rels)
        
        # Getting the group algebra over GF(2)
        F2G = group_algebra(GF(2), Q)
        
        # Grab the generators so we can use them later.
        q_gens = gens(Q)
        return Q, F2G, q_gens[1], q_gens[2]
    end

    @testset "Simple ZSZ Code Construction" begin
        # Setting a seed so this doesn't flake out on us in CI.
        Random.seed!(42)

        # I'm using the example Z_5 x| Z_4 with q=2.
        # Checked that 2^4 = 16 which is 1 mod 5, so we're good.
        # This gives us a group of size 20.
        l, m, q = 5, 4, 2
        
        Q, F2G, x, y = zsz_group_rings(l, m, q)
        id = one(Q)

        # The paper mentions searching for polynomials.
        # I'm gonna do a probabilistic search here for 3-term polynomials.
        # A = 1 + ... + ...
        # B = 1 + ... + ...
        
        found_good_code = false
        
        # Let's give it 20 shots to find something decent.
        for _ in 1:20
            # Helper to just grab unique random elements like x^i * y^j
            function rand_poly_elts()
                elts = [id]
                while length(elts) < 3
                    # Picking random powers for x and y
                    i = rand(0:l-1)
                    j = rand(0:m-1)
                    elt = x^i * y^j
                    if !(elt in elts)
                        push!(elts, elt)
                    end
                end
                return elts
            end

            # Generating my random polynomials A and B
            A_elts = rand_poly_elts()
            B_elts = rand_poly_elts()

            # Constructing the 2BGA code with them
            c = twobga_from_fp_group(A_elts, B_elts, F2G)
            
            # Now let's check if it actually worked.
            # n should be 2 * size of group (2 * 20 = 40).
            n = code_n(c)
            k = code_k(c)
            
            # I want a code that actually encodes something (k > 0).
            if n == 2*l*m && k > 0
                @test n == 40
                @test k >= 2 
                found_good_code = true
                
                # If we had a fast solver, I'd check distance > 1 too, but skipping for speed.
                break
            end
        end
        
        if !found_good_code
            @warn "Didn't find a valid ZSZ code in my attempts. Maybe verify the seed?"
        end
        @test found_good_code
    end
    
    @testset "Invalid ZSZ Parameters" begin
        # Just making sure I catch bad parameters.
        # 2^3 = 8 which is 3 mod 5, not 1. Should fail.
        @test_throws ArgumentError zsz_group_rings(5, 3, 2)
    end
end
