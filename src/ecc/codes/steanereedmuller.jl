"""The family of non-CSS Steane-Reed-Muller codes, as discovered by Steane in his 1999 paper [steane1999quantum](@cite).

You might be interested in consulting [zhang1997quantum](@cite), [quan2018fault](@cite) and [campbell2012magic](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_reed_muller)
"""
struct SteaneReedMuller <: AbstractECC
    t::Int  
    r::Int  
    function SteaneReedMuller(t, r)
        if t == 1 && r == 2 
            throw(ArgumentError("Invalid parameters when attempting to construct a Steane-Reed-Muller code. t = 1, r = 2 is not quantum code as it does not satisfy commutativity requirements."))
        elseif t < 0 || t > 7 || r < 1 || r > 7
            throw(ArgumentError("Invalid parameters when attempting to construct a Steane-Reed-Muller code. We need  0<r<7 and 0≤t<7 in order to obtain a valid code and to remain tractable."))
        else
            new(t, r)
        end
    end
end

function _codesize_qrm(t, r)
    return 2^r - sum(binomial.(r, 0:t)) 
end 

function _nplusk_qrm(t, r)
    return (2^r + 2^r) - sum(binomial.(r, 0:t)) -  sum(binomial.(r, 0:t - 1)) 
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
    if Dx_rows - r == 0
        J = parity_checks(ReedMuller(t + 1, r))
        J = _steane_convention_qrm(J)
        Dx = J[end-r+1:end, :]
    elseif _nplusk_qrm(t, r) - 2^r == 0
        G1 = G1[1:_codesize_qrm(t, r), :]
        K = parity_checks(ReedMuller(t, r))
        K = _steane_convention_qrm(K)
        Dx = K[end-Dx_rows+1:end, :]
    else
        Dx = G1[end-Dx_rows+1:end, :]
    end
    #Dz
    m = _codesize_qrm(t - 1, r) - _codesize_qrm(t, r)
    Dz = zeros(Int64, size(Dx, 1), size(Dx, 2))
    for i in 1:m - 1
        Dz[i, :] = Dx[i + 1, :]
    end
    e = (_codesize_qrm(t, r) - Dx_rows - t) == 0 || (_codesize_qrm(t, r) - Dx_rows - t) < 0 ? t : _codesize_qrm(t, r) - Dx_rows - t
    Dz[m, :] = circshift(Dx[1, :], -e) 
    
    Hx = vcat(G1, zeros(Int64, size(G1, 1), size(G1, 2)), Dx)
    Hz = vcat(zeros(Int64, size(G1, 1), size(G1, 2)), G1, Dz)
    extended_Hx = Matrix{Bool}(Hx)
    extended_Hz = Matrix{Bool}(Hz)
    Stabilizer(extended_Hx, extended_Hz)
end