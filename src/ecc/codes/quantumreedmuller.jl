"""The family of quantum Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller)
"""
struct QuantumReedMuller <: AbstractECC
    t::Int  
    r::Int  
    function QuantumReedMuller(t, r)
        if t == 1 && r == 2 
            throw(ArgumentError("Invalid parameters: t = 1, r = 2 is not quantum code as it does not satisfy Hx.Hz' + Hz.Hx' = 0 [steane1999quantum](@cite)."))
        elseif t < 0 || t > 10 || r < 1 || r > 10
            throw(ArgumentError("Invalid parameters: r must be positive and < 10 and t >= 0 and < 10 in order to obtain a valid code and to remain tractable."))
        elseif 2^t + 2^(t + 1) > 2^r
            throw(ArgumentError("Invalid parameters: The minimum Hamming distance (2^t + 2^(t + 1)) must be < 2^r to obtain a valid code and to remain tractable."))
        else
        new(t, r)
        end
    end
end

function iscss(::Type{QuantumReedMuller})
    return true
end

function generate_Dx(dx_rows, row_len, hdist)
    combs = []
    rows = 0

    for j in combinations(1:row_len, hdist)
        comb = zeros(Int, row_len)
        if rows >= dx_rows
             break
        end
        for k in j
            comb[k] = 1
        end
        push!(combs, comb)
        rows += 1
    end
    r = length(combs)
    c = length(combs[1])
    H = reshape(vcat(combs...), c, r)'
    return H
end

function trows(t, r)
    return (2^r + 2^r) - sum(binomial.(r, 0:t)) -  sum(binomial.(r, 0:t - 1)) 
end 

function k0(t, r)
    return 2^r - binomial(r, t) -  sum(binomial.(r, 0:t - 1)) 
end 

function parity_checks(c::QuantumReedMuller)
    r = c.r
    t = c.t
    
    G1 = parity_checks(ReedMuller(t, r))
    pad_zeros = zeros(Int64, size(G1, 1), size(G1, 2))

    total_rows = trows(t, r)
    if (total_rows - 2^r) == 0
        Dx_rows = r + t
        k = k0(t, r)
        G1 = G1[1:k, :]
        pad_zeros = zeros(Int64, size(G1, 1), size(G1, 2))     
    else
        Dx_rows = abs(total_rows - 2^r)
    end
   
    Dx_cols = size(G1, 2)
    Dx = generate_Dx(Dx_rows, Dx_cols, 2^t)
   
    m = sum(binomial.(r, 0:t)) - sum(binomial.(r, 0:t - 1))
    Dz = zeros(Int64, size(Dx, 1), size(Dx, 2))
    for i in 1:m - 1
        Dz[i, :] = Dx[i + 1, :]
    end
    Dz[m, :] = circshift(Dx[1, :], -m) 
    
    Hx = vcat(G1, pad_zeros, Dx)
    Hz = vcat(pad_zeros, G1, Dz)  

    extended_Hx = Matrix{Bool}(Hx)
    extended_Hz = Matrix{Bool}(Hz)

    num_rows = size(Hx, 1)
    fill_array = fill(UInt8(0), num_rows)
    Stabilizer(fill_array, extended_Hx, extended_Hz)
end