GPU_DEPENDENCIES = ["CUDA"]

"""
copies the memory content of the object to GPU
"""
function to_cpu end

"""
copies the memory content of the object to CPU
"""
function to_gpu end
