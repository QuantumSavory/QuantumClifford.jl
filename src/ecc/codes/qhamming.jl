"""This class implements the Gottesman codes [[2^j, 2^j-j-2,3]], also known as quantum Hamming codes, as described in Gottesman's 1997 paper on stabilizer codes [gottesman1997stabilizer](@cite)."""

struct QHamming <: AbstractECC
    j::Int
    function QHamming(j)
        (j >= 3 && j < 21) || error("In `QHamming(j)`, `j` must be ≥  3 in order to obtain a valid code and `j` must be < 21 to remain tractable")
        new(j)
    end
end

code_n(c::QHamming)= 2^c.j

function parity_checks(c::QHamming)
    rows = c.j + 2
    cols = 2^c.j
    
    Hx = falses(rows,cols)
    Hz = falses(rows,cols)
    
    Hx[1, :] .= true
    
    Hx[2, :] .= false

    for i in 3:(rows-1)
        for a in 1:cols
            Hx[i, a] = (a==0) || (a % 2 ==0)
        end
    end
   
    for a in 1:cols
        Hx[rows, a] = (a % 4 == 1) || (a % 4 == 2) ? 0 : 1
    end
   
    Hz[1, :] .= false
    Hz[2, :] .= true
    
    for i in 3:rows
        period = 2^(rows - i)
        for a in 1:cols
            Hz[i,a]=div(a-1,period) % 2 == 1
        end
    end
    extended_Hx = Matrix{Bool}(Hz)
    extended_Hz = Matrix{Bool}(Hx)
    
    num_rows = size(Hx, 1)
   
    fill_array = fill(UInt8(0), num_rows)
    Stabilizer(fill_array, extended_Hz, extended_Hx)
end