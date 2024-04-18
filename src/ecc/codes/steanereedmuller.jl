"""The family of non-CSS Steane-Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

The construction method is effective for k > 0, and generators G  = (G_x | G_z) for k < 0 serve as stabilizers (parity check matrices) H = (H_x | H_z) for k > 0 codes.

You might be interested in consulting [zhang1997quantum](@cite), [quan2018fault](@cite), and [campbell2012magic](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller).
"""
struct SteaneReedMuller <: AbstractECC
    t::Int  
    r::Int  
    function SteaneReedMuller(t, r)
        if t < 0 || t > 7 || r < 1 || r > 7
            throw(ArgumentError("Invalid parameters when attempting to construct a Steane-Reed-Muller code. We need  0<r≤7 and 0≤t≤7 in order to obtain a valid code and to remain tractable."))
        elseif (2^r - binomial(r, t) - 2*sum(binomial.(r, 0:t - 1)))  < 0 || (2^r - binomial(r, t) - 2*sum(binomial.(r, 0:t - 1)))  == 0
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

function _steane_convention_qrm(mat)
    inv_mat = copy(mat)
    for r in 1:size(mat, 1)
        inv_mat[r, :] = reverse(mat[r, :])
    end
    return inv_mat
end

function parity_checks(c::SteaneReedMuller)
    r = c.r
    t = c.t
    #G1
    G1 = parity_checks(ReedMuller(t, r))
    G1 = _steane_convention_qrm(G1)
    #Dx
    Dx_rows = _codesize_qrm(t - 1, r) - _codesize_qrm(t, r)
    # Gottesman Codes, @ t = 1, r = 3 ... 9
    if Dx_rows - r == 0
        J = parity_checks(ReedMuller(t, r))
        J = _steane_convention_qrm(J)
        Hx =  J[1:1, :]
        Hx = vcat(Hx, zeros(Int64, 1, size(J, 2)))
        Hx = vcat(Hx, J[2:end,:])
        sa = div(length(J[2, :]), 4)
        sr = circshift(J[2, :], -sa)
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
    Dz[m, :] = circshift(Dx[1, :], -t^2) 
    pad_zeros = zeros(Int64, size(G1, 1), size(G1, 2))
    Hx = zeros(Int64, size(G1, 1) + size(pad_zeros, 1) + size(Dx, 1), size(G1, 2))
    Hz = zeros(Int64, size(G1, 1) + size(pad_zeros, 1) + size(Dz, 1), size(G1, 2)) 
    Hx = vcat(G1, pad_zeros, Dx)
    Hz = vcat(pad_zeros, G1, Dz) 
    bool_Hx = Matrix{Bool}(Hx)
    bool_Hz = Matrix{Bool}(Hz)
    Stabilizer(bool_Hx, bool_Hz)
end
