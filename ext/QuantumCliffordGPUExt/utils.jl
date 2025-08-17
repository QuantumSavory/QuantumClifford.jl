const THROW_INVALID_CUDA_ARG = 
"First argument to @run_cuda should be a function call"
macro run_cuda(call, ndrange)
    # destructure the kernel call
    Meta.isexpr(call, :call) || throw(ArgumentError(THROW_INVALID_CUDA_ARG))
    f = call.args[1]
    args = call.args[2:end]
    args = [esc(x) for x in args]
    quote
        local kernel = @cuda launch=false $f($(args...))
        local config = CUDA.launch_configuration(kernel.fun)
        local threads = min($(esc(ndrange)), config.threads)
        local blocks = cld($(esc(ndrange)), threads)
        kernel($(args...); threads=threads, blocks=blocks)
    end
end
