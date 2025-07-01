include("test_KA_mul_leftright.jl")

@inline function test_platform(AT, synchronize)
	@testset "mul_leftright" begin
		test_KA_mul_leftright(AT, synchronize)
	end
end
