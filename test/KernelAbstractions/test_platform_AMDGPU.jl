@testitem "AMDGPU" tags = [:amdgpu] begin

	using AMDGPU: ROCArray
	AT = ROCArray
	import KernelAbstractions as KA
	backend = KA.get_backend(AT([0]))
	synchronize() = KA.synchronize(backend)

	@testset "mul_leftright" begin
		include("test_KA_mul_leftright.jl")
		test_KA_mul_leftright(AT, synchronize)
	end

end
