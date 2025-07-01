@testitem "OpenCL" begin

	include("test_platform.jl")

	import pocl_jll
	using OpenCL: CLArray, cl.devices, cl.platforms, cl.finish, cl.queue
	AT = CLArray

	can_run = any(length(devices(platform)) > 0 for platform in platforms())

	@testset "Device availability" begin
		@test can_run
	end

	if can_run
		synchronize() = finish(queue())
		test_platform(AT, synchronize)
	end

end
