include("benchmark_KA_mul_leftright.jl")

@inline function benchmark_platform(
	AT, synchronize;
	platform_name = string(AT)
	)

	benchmark_KA_mul_leftright(
		AT, synchronize;
		phases = Val(true), platform_name = platform_name
		)
	benchmark_KA_mul_leftright(
		AT, synchronize;
		phases = Val(false), platform_name = platform_name
		)

end
