abstract type ClassicalCode end

struct ReedMuller <: ClassicalCode
    r::Int
    m::Int

    function ReedMuller(r, m)
        if r < 0 || m < 1
            throw(ArgumentError("Invalid parameters: r must be non-negative and m must be positive"))
        end
        new(r, m)
    end
end

function construct(m, i)
    return repeat([fill(1, 2^(m - i - 1)); fill(0, 2^(m - i - 1))], outer = 2^i)
end

function vmult(vecs...)
    return [reduce(*, a, init=1) for a in zip(vecs...)]
end

function parity_checks(c::ReedMuller)
    r=c.r
    m=c.m
    rx = [construct(m, i) for i in 0:m - 1]
    row_matrices = [reduce(vmult, [rx[i + 1] for i in S], init = ones(Int, 2^m)) 
                     for s in 0:r for S in combinations(0:m - 1, s)]
    rows = length(row_matrices)
    cols = length(row_matrices[1])
    H = reshape(vcat(row_matrices...), cols, rows)'
end