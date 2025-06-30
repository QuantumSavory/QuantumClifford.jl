using AMDGPU: ROCArray
AT = ROCArray
import KernelAbstractions as KA
backend = KA.get_backend(AT([0]))
synchronize() = KA.synchronize(backend)

include("benchmark_KA_mul_leftright.jl")
benchmark_mul_leftright(
	AT, synchronize;
	phases = Val(true), platform_name = "AMDGPU"
	)
benchmark_mul_leftright(
	AT, synchronize;
	phases = Val(false), platform_name = "AMDGPU"
	)
