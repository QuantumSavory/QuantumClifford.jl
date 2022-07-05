function test_hash()
    @testset "Hashing" begin
        @test hash(P"X") == hash(P"X")
        @test hash(S"X") == hash(S"X")
        @test hash(C"X Z") == hash(C"X Z")
    end
end

test_hash()
