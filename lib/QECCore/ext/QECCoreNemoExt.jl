module QECCoreNemoExt

using QECCore
using DocStringExtensions

import Nemo
import Nemo: GF, matrix, rank, transpose

import QECCore: TillichZemor, code_k, _create_circulant_matrix, _create_matrix_M_deterministic

function QECCore.code_k(c::TillichZemor)
    C = _create_circulant_matrix(c.m)
    M = _create_matrix_M_deterministic(c.m, c.n, c.r)
    H = hcat(C, M)
    H_gf2 = matrix(GF(2), H)
    Ht_gf2 = transpose(H_gf2)
    k = c.n - rank(H_gf2)
    kT = c.m - rank(Ht_gf2)  # c.m == size(H, 1)
    k_q = k^2 + kT^2
    return k_q
end

end # module
