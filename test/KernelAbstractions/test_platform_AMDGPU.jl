@testitem "AMDGPU" tags = [:amdgpu] begin

	include("test_platform.jl")

	using AMDGPU: ROCArray, synchronize, devices
	AT = ROCArray

	can_run = length(devices()) > 0

	@testset "Device availability" begin
		@test can_run
	end

	if can_run
		test_platform(AT, synchronize)
	end

end
