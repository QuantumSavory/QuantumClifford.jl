module QECCore

# interfaces
export distance, parity_matrix, code_n, code_s, code_k, parity_matrix_x, parity_matrix_z, rate
export AbstractECC, AbstractQECC, AbstractCECC, AbstractCSSCode

include("interface.jl")

function __init__()
    Base.Experimental.register_error_hint(_hint_MethodError, MethodError)
end

function _hint_MethodError(io::IO, exc::MethodError)
    if exc.f == parity_matrix_x || exc.f == parity_matrix_z
        print(io, "\n╰─ Hint: This error usually occurs when trying to call `parity_matrix_x` or `parity_matrix_z` on a non-CSS quantum code.")
    end
end
end
