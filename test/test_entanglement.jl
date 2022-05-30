function test_clipping()
    @testset "Clip gauge of stabilizer states" begin
        for n in test_sizes
            s = random_stabilizer(n)
            s_clipped = copy(s)
            clip!(s_clipped)
            @test logdot(s, s_clipped)==0
            bigram = get_bigram(s_clipped; do_clip=false)
            rows, columns = size(stabilizerview(s_clipped))
            for j in 1:columns
                @test count(==(j), bigram)==2
                #TODO: we do not whether test endpoints are different Pauli operators
            end
        end
    end
end

test_clipping()