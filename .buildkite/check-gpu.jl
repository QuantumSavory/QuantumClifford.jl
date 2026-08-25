backend = only(ARGS)

if backend == "cuda"
    import CUDA
    CUDA.functional() || error("CUDA is not functional on this runner")
elseif backend == "rocm"
    import AMDGPU
    isempty(AMDGPU.devices()) && error("No ROCm device is available on this runner")
elseif backend == "opencl"
    import OpenCL, pocl_jll
    any(!isempty(OpenCL.cl.devices(platform)) for platform in OpenCL.cl.platforms()) ||
        error("No OpenCL device is available on this runner")
else
    error("Unknown GPU backend: $backend")
end
