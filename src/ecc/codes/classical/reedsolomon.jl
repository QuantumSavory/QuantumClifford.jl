
struct ReedSolomon <: ClassicalCode
    n::Int
    k::Int
    j::Int

    function ReedSolomon(n, k, j)
        if n < 0 || k < 0 || j < 0 || n >= 500
            throw(ArgumentError("Invalid parameters: n, k, j must be non-negative and n >= 500 in order to obtain a valid code and to remain tractable"))
        end 
            new(n, k, j)
    end
end

function generator_matrix(rs::ReedSolomon)
    r = ceil(Int, log2(rs.n + 1))
    t = div(rs.n - rs.k, 2)
    GF2ͬ, a = finite_field(2, r, "a")
    P, x = GF2ͬ[:x]
    polynomial = x - a^rs.j
    for i in 1:(2*t - 1)
       polynomial *= (x - a^(rs.j + i))
    end
    return polynomial
end
