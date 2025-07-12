@testitem "ROCm" tags = [:rocm] begin

	include("test_platform.jl")

	using AMDGPU: ROCArray, synchronize, devices
	const AT = ROCArray

	const can_run = length(devices()) > 0

	@testset "Device availability" begin
		@test can_run
	end

	if can_run
		test_platform(AT, synchronize)
	end

end
