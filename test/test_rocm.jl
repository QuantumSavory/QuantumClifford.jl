@testitem "ROCM Tests" tags=[:rocm] begin
    using AMDGPU
    arr = rand(1000)
    @test ROCArray(arr)+ROCArray(arr) == ROCArray(arr+arr)
end
