"""The family of quantum Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller)
"""
struct QuantumReedMuller <: AbstractECC
    t::Int  
    r::Int  
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
        rows +=1
    end
    r = length(combs)
    c = length(combs[1])
    H = reshape(vcat(combs...), c, r)'
    return H
end

function trows(t, r)
    return 2^r + 2^r - sum(binomial.(r, 0:t)) -  sum(binomial.(r, 0:t - 1)) 
end 

function parity_checks(c::QuantumReedMuller)
    r = c.r
    t = c.t
    
    G1 = parity_checks(ReedMuller(t, r))
    pad_zeros = zeros(Int64, size(G1, 1), size(G1, 2))
     
    total_rows = trows(t, r)
    if t == 2 && r == 1
        Dx_rows = 0     
    else
        Dx_rows = abs(total_rows - 2^r)
    end

    Dx_cols = size(G1, 2)
    Dx = generate_Dx(Dx_rows, Dx_cols, t)
    
    m = sum(binomial.(r, 0:t)) - sum(binomial.(r, 0:t - 1))
    Dz = Dx
    for i in 1:m - 1
        Dz[i, :] .= Dx[i + 1, :]
    end
    Dz[m, :] .= reverse(circshift(reverse(Dx[1, :]), 1))

    Hx = vcat(G1, pad_zeros, Dx)    
    Hz = vcat(pad_zeros, G1, Dz)    

    extended_Hx = Matrix{Bool}(Hx)
    extended_Hz = Matrix{Bool}(Hz)

    num_rows = size(Hx, 1)
    fill_array = fill(UInt8(0), num_rows)
    Stabilizer(fill_array, extended_Hx, extended_Hz)
end