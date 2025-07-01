@testitem "CUDA" tags = [:cuda] begin

	include("test_platform.jl")

	using CUDA: CuArray, synchronize, devices
	AT = CuArray

	can_run = length(devices()) > 0

	@testset "Device availability" begin
		@test can_run
	end

	if can_run
		test_platform(AT, synchronize)
	end

end
