GPU_DEPENDENCIES = ["CUDA"]

# todo this is problematic because the error message is not specific
# todo make a better error stack trace?
_throw_error_gpu_dependencies(T::Type) = error("either this function is not supported on type $T or dependencies are not imported.\ndependencies: $GPU_DEPENDENCIES")

to_cpu(gpu_object::T) where {T} = _throw_error_gpu_dependencies(T)
to_gpu(cpu_object::T) where {T} = _throw_error_gpu_dependencies(T)
