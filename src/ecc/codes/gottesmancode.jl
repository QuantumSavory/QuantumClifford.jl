"""This class implements the family of [[2^j, 2^j - j - 2, 3]] Gottesman codes, also known as [[2^j, 2^j - j - 2, 3]] quantum Hamming codes, as described in Gottesman's 1997 PhD thesis on 'Stabilizer Codes and Quantum Error Correction' [Gottesman1997](@cite).

Key details for the [[8, 3, 3]] code when j = 3 using Stabilizer Formalism:
- [6506105](@cite) utilizes an explicit construction method for the stabilizer generators. The remaining m generators, apart from Mx and Mz, are constructed using the check matrix [Hm|AmHm], where Hm = [c0, c1, ..., c^2m − 1]. Each column ck (for k = 0, 1, ..., 2^m − 1) represents a binary vector corresponding to the integer k. Also, Am refers to any invertible m × m matrix devoid of fixed points such that Am​.s ≠ 0 and Am.s ≠ s for all s ∈ F₂ᵐ.
- [Chao2017QuantumEC](@cite) introduces alternative stabilizer generators for the [[8, 3, 3]] code. It identifies permutations effective for all stabilizer generators except for X^⊗8 and Z^⊗8. However, these can be replaced with XXYZIYZI and ZZIXYIXY, respectively. The stabilizer generators are permuted to achieve desired properties, reflecting an alternative approach to constructing the stabilizer code.

Notes:
- This implementation adopts the Gottesman Notation for the [[8, 3, 3]] stabilizer code as depicted in Table 3.3 on Page 22 on the Gottesman's PhD thesis. Additionally, for j = 4, our [[16, 10, 3]] stabilizer code stemming from QHamming(4) can be cross-verified from Table 8.1 on Page 91 as well [Gottesman1997](@cite).
- The differences between [6506105](@cite), [Chao2017QuantumEC](@cite), and [Gottesman1997](@cite) primarily lie in the choice of stabilizer generators and their permutations for the [[8, 3, 3]] stabilizer code.
- The discrepancy in stabilizer generator representations underscores the flexibility within the stabilizer formalism, allowing for various valid choices and permutations, albeit with different notations and implementation details.
"""

struct Gottesman <: AbstractECC
    j::Int
    function Gottesman(j)
        (j >= 3 && j < 21) || error("In `Gottesman(j)`, `j` must be ≥  3 in order to obtain a valid code and `j` must be < 21 to remain tractable")
        new(j)
    end
end

code_n(c::Gottesman) = 2^c.j

function parity_checks(c::Gottesman)
    rows = c.j + 2
    cols = 2^c.j
    
    Hx = falses(rows,cols)
    Hz = falses(rows,cols)
    
    Hx[1, :] .= true
    Hx[2, :] .= false

    if c.j == 3 
        for col in 1:cols
            Hx[3, col] = (col % 8 == 1 || col % 8 == 3 || col % 8 == 6) ? 0 : 1
        end
        Hx[3, cols] = Hx[3, cols] == 0 ? 1 : 0
        for col in 1:cols
            Hx[4, col] = (col % 4 == 1) || (col % 4 == 3) ? 0 : 1
        end
        for a in 1:cols
            Hx[rows, a] = ((a % 4 == 0) || (a % 4 == 1) ? 0 : 1) ⊻ ((a % 8 == 5) || (a % 8 == 6))
        end
        Hx[end, [end-1, end]] .= [0, 1]

    else
        for a in 1:cols
            Hx[3, a] = (a == 0) || (a % 2 == 0)
        end
        for row in 4:rows - 1
            for col in 1:cols
                k = row - 3
                m = 2^(c.j - k)
                n = 2^(c.j - k)
                if (col - 1) % (m + n) < m
                    if col % 2 == 0
                        Hx[row, col] = 1
                    else
                        Hx[row, col] = 0
                    end
                else
                    if col % 2 == 0
                        Hx[row, col] = 0
                    else
                        Hx[row, col] = 1
                    end
                end
            end
        end
    
        for a in 1:cols
            Hx[rows, a] = (a % 4 == 0) || (a % 4 == 1) ? 0 : 1
        end
    end

    Hz[1, :] .= false
    Hz[2, :] .= true
    
    for i in 3:rows
        period = 2^(rows - i)
        for a in 1:cols
            Hz[i, a] = div(a - 1, period) % 2 == 1
        end
    end
    extended_Hx = Matrix{Bool}(Hz)
    extended_Hz = Matrix{Bool}(Hx)
    
    num_rows = size(Hx, 1)
   
    fill_array = fill(UInt8(0), num_rows)
    Stabilizer(fill_array, extended_Hz, extended_Hx)
end