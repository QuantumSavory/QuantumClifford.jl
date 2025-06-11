module QECCore

# interfaces
export distance, parity_matrix, code_n, code_s, code_k, parity_matrix_x, parity_matrix_z, rate
export AbstractECC, AbstractQECC, AbstractCECC, AbstractCSSCode, AbstractDistanceAlg

include("interface.jl")

function __init__()
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
            if exc.f == parity_matrix_x || exc.f == parity_matrix_z
                print(io, "\n Hint: This error usually occurs when trying to call `parity_matrix_x` or `parity_matrix_z` on a non-CSS quantum code.")
            end
        end
    end
end

end
