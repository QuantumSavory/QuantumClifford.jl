GPU_DEPENDENCIES = ["CUDA"]

# todo this is problematic because the error message is not specific
_throw_error_gpu_dependencies(T::Type) = error("either or this function is not supported on type $T on dependencies are not imported.\ndependencies: $GPU_DEPENDENCIES")

to_cpu(gpu_object::T) where {T} = _throw_error_gpu_dependencies(T)
to_gpu(cpu_object::T) where {T} = _throw_error_gpu_dependencies(T)
