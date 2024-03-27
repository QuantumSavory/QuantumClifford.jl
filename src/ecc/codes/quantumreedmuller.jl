"""The family of quantum Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

You might be interested in consulting [zhang1997quantum](@cite), [quan2018fault](@cite) and [campbell2012magic](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller)
"""
struct QuantumReedMuller <: AbstractECC
    t::Int  
    r::Int  
    function QuantumReedMuller(t, r)
        if t == 1 && r == 2 
            throw(ArgumentError("Invalid parameters when attempting to construct a quantum Reed-Muller code. t = 1, r = 2 is not quantum code as it does not satisfy commutativity requirements."))
        elseif t < 0 || t > 7 || r < 1 || r > 7
            throw(ArgumentError("Invalid parameters: r must be positive and < 7 and t >= 0 and < 7 in order to obtain a valid code and to remain tractable."))
        else
            new(t, r)
        end
    end
end

function codesize(t, r)
    return 2^r - sum(binomial.(r, 0:t)) 
end 

function nplusk(t, r)
    return (2^r + 2^r) - sum(binomial.(r, 0:t)) -  sum(binomial.(r, 0:t - 1)) 
end 

function parity_checks(c::QuantumReedMuller)
    r = c.r
    t = c.t
    
    G1 = parity_checks(ReedMuller(t, r))
    
    Dx_rows = codesize(t - 1, r) - codesize(t, r)
    if Dx_rows - r == 0
        J = parity_checks(ReedMuller(t + 1, r))
        Dx = J[end-r+1:end, :]
    elseif nplusk(t, r) - 2^r == 0
        G1 = G1[1:codesize(t, r), :]
        K = parity_checks(ReedMuller(t, r))
        Dx = K[end-Dx_rows+1:end, :]
    else
        Dx = G1[end-Dx_rows+1:end, :]
    end
    
    m = codesize(t - 1, r) - codesize(t, r)
    Dz = zeros(Int64, size(Dx, 1), size(Dx, 2))
    for i in 1:m - 1
        Dz[i, :] = Dx[i + 1, :]
    end
    e = (codesize(t, r) - Dx_rows - t) == 0 || (codesize(t, r) - Dx_rows - t) < 0 ? t : codesize(t, r) - Dx_rows - t
    
    Dz[m, :] = circshift(Dx[1, :], e) 
   
    Hx = vcat(G1, zeros(Int64, size(G1, 1), size(G1, 2)), Dx)
    Hz = vcat(zeros(Int64, size(G1, 1), size(G1, 2)), G1, Dz)  

    extended_Hx = Matrix{Bool}(Hx)
    extended_Hz = Matrix{Bool}(Hz)

    num_rows = size(Hx, 1)
    fill_array = fill(UInt8(0), num_rows)
    Stabilizer(fill_array, extended_Hx, extended_Hz)
end