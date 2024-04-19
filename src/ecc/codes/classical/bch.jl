"""The family of Bose–Chaudhuri–Hocquenghem (BCH) codes, as discovered in 1959 by Alexis Hocquenghem [hocquenghem1959codes](@cite), and independently in 1960 by Raj Chandra Bose and D.K. Ray-Chaudhuri [bose1960class](@cite).

You might be interested in consulting [bose1960further](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/q-ary_bch)
"""
struct BCH <: ClassicalCode
    n::Int 
    t::Int #Error correction capability; t bits can be corrected

    function BCH(n, t)
        if n < 6 || n > 500 || t < 0 || t > 2^(ceil(Int, log2(n + 1)) - 1)
            throw(ArgumentError("Invalid parameters: 'n' and 't' must be positive, and 'r' must be >= to 3. Additionally, 'n' is >= to 7 since n = 2ͬ - 1 and 't' < 2^(r - 1), to obtain a valid code and to tractable."))

        end
            new(n, t)
    end
end

function generator_polynomial(rs::BCH)
    r = ceil(Int, log2(rs.n + 1))
    GF2ͬ, a = finite_field(2, r, "a")
    GF2x, x = GF(2)["x"]
    minimal_poly = FqPolyRingElem[]
    for i in 1:(2*rs.t - 1)
        if i % 2 != 0
             minimal_poly = [minimal_poly; minpoly(GF2x, a^i)]
        end 
    end
    gx = lcm(minimal_poly)
    return gx
end

function parity_checks(rs::BCH)
    r = ceil(Int, log2(rs.n + 1))
    GF2ͬ, a = finite_field(2, r, "a")
    HField = Matrix{FqFieldElem}(undef, rs.t, rs.n)
    for i in 1:rs.t
        for j in 1:rs.n
            base = 2*i - 1  
            HField[i, j] = base^(j-1) * a^(j-1)
        end
    end
    H = Matrix{Bool}(undef, r*rs.t, rs.n)
    for i in 1:rs.t
        row_start = (i - 1) * r + 1
        row_end = row_start + r - 1
        for j in 1:rs.n
            t_tuple = Bool[]
            for k in 0:r - 1
                t_tuple = [t_tuple; !is_zero(coeff(HField[i, j], k))]
            end 
            H[row_start:row_end, :] .=  vec(t_tuple')
        end
    end 
    return H
end
