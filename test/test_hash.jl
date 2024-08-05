@testitem "Hashing" begin
    @test hash(P"X") == hash(P"X")
    @test hash(S"X") == hash(S"X")
    @test hash(C"X Z") == hash(C"X Z")
    r = random_destabilizer(2)
    @test hash(r) == hash(copy(r))
end
