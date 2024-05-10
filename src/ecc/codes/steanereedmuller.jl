"""The family of non-CSS Steane-Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

The generators G = (G_x | G_z) of Steane-Reed-Muller codes are constructed using the classical Reed-Muller codes where RS(t, r) primitive generates all the evaluation vectors of polynomials with t variables and degree no larger than r.

The construction method is effective for code size k > 0, and generators G for k < 0 serve as stabilizers (parity check matrices) H = (H_x | H_z) for k > 0 codes. This is because the stabilizer (H) of a code [[n, k, d]] is equal to the generator matrix (G) of its dual code [[n, -k, d']]. The relationship between the original code's minimum distance (d) and the dual code's minimum distance (d') is determined by the parity (odd or even) of parameter r, as detailed in Table 1 [steane1999quantum](@cite). This equivalence between H and the dual's G allows these codes satisfy "self-dual" condition. 

You might be interested in consulting [zhang1997quantum](@cite), [quan2018fault](@cite), and [campbell2012magic](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller).
"""
struct SteaneReedMuller <: AbstractECC
    t::Int  
    r::Int  
    function SteaneReedMuller(t, r)
        if t < 0 || t > 7 || r < 1 || r > 7
            throw(ArgumentError("Invalid parameters when attempting to construct a Steane-Reed-Muller code. We need  0<r≤7 and 0≤t≤7 in order to obtain a valid code and to remain tractable."))
        elseif (2^r - binomial(r, t) - 2*sum(binomial.(r, 0:t - 1)))  <= 0
            throw(ArgumentError("Invalid parameters when attempting to construct a Steane-Reed-Muller code. The method to construct Steane-Reed-Muller code fails when k < 0."))  
        else
            new(t, r)
        end
    end
end

#Equation 10 [steane1999quantum](@cite)
function _codesize_qrm(t, r)
    return 2^r - sum(binomial.(r, 0:t)) 
end 

function parity_checks(c::SteaneReedMuller)
    r = c.r
    t = c.t
    #G1
    G1 = parity_checks(ReedMuller(t, r))
    #Dx
    Dx_rows = _codesize_qrm(t - 1, r) - _codesize_qrm(t, r)
    # Gottesman Codes, @ t = 1, r = 3 ... 9
    if Dx_rows - r == 0
        J = parity_checks(ReedMuller(t, r))
        Hx =  J[1:1, :]
        Hx = vcat(Hx, zeros(Int64, 1, size(J, 2)))
        Hx = vcat(Hx, J[2:end,:])
        sa = div(length(J[2, :]), 4)
        sr = circshift(J[2, :], sa)
        sr = reshape(sr, 1, length(sr))  
        Hz = zeros(Int64, size(J, 1), size(J, 2))         
        Hz[2, :] = J[1:1, :]
        Hz[3:end, :] = J[3:end, :]
        Hz = vcat(Hz, sr)
        gHx = Matrix{Bool}(Hx)
        gHz = Matrix{Bool}(Hz)
        return Stabilizer(gHx, gHz)
    else
        Dx = zeros(Int64, Dx_rows, size(G1, 2))
        Dx = G1[end-Dx_rows+1:end, :]
    end
    #Dz
    m = _codesize_qrm(t - 1, r) - _codesize_qrm(t, r)
    Dz = zeros(Int64, size(Dx, 1), size(Dx, 2))
    for i in 1:m - 1
        Dz[i, :] = Dx[i + 1, :]
    end
    Dz[m, :] = circshift(Dx[1, :], t^2) 
    pad_zeros = zeros(Int64, size(G1, 1), size(G1, 2))
    Hx = vcat(G1, pad_zeros, Dx)
    Hz = vcat(pad_zeros, G1, Dz) 
    bool_Hx = Matrix{Bool}(Hx)
    bool_Hz = Matrix{Bool}(Hz)
    Stabilizer(bool_Hx, bool_Hz)
end
