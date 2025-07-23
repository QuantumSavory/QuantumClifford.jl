module QECCoreNemoExt

using QECCore
using DocStringExtensions

import Nemo
import Nemo: GF, matrix, rank, transpose

import QECCore: code_k, parity_matrix_x, parity_matrix_z

function QECCore.code_k(c::AbstractQECC)
    n = code_n(c)
    F₂ = GF(2)
    rank_Hx = rank(matrix(F₂, parity_matrix_x(c)))
    rank_Hz = rank(matrix(F₂, parity_matrix_z(c)))
    return n - rank_Hx - rank_Hz
end

end # module
