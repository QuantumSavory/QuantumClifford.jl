@testitem "CUDA Tests" tags=[:cuda] begin
    using CUDA
    arr = rand(1000)
    @test CuArray(arr)+CuArray(arr) == CuArray(arr+arr)
end