"""The family of `[[2ʲ, 2ʲ - j - 2, 3]]` Gottesman codes, also known as quantum Hamming codes, as described in [Gottesman's 1997 PhD thesis](@cite Gottesman1997).

You might be interested in consulting [yu2013all](@cite) and [chao2018quantum](@cite) as well.

The ECC Zoo has an [entry for this family](https://errorcorrectionzoo.org/c/quantum_hamming)
"""
struct Gottesman <: AbstractECC
    j::Int
    function Gottesman(j)
        (j >= 3 && j < 21) || error("In `Gottesman(j)`, `j` must be ≥  3 in order to obtain a valid code and < 21 to remain tractable")
        new(j)
    end
end

code_n(c::Gottesman) = 2^c.j

function parity_checks(c::Gottesman)
    rows = c.j + 2
    cols = 2^c.j

    Hx = falses(rows, cols)
    Hz = falses(rows, cols)

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

        if rows % 2 == 1
            for a in 1:cols
                Hx[rows, a] = ((a % 4 == 0) || (a % 4 == 1) ? 0 : 1) ⊻ ((a % 8 == 5) || (a % 8 == 6) || (a % 8 == 7))
            end
               for a in 1:div(cols, 8)
                   ci = 8 * a
                   if ci <= cols
                       Hx[rows, ci] = 1 -Hx[rows, ci]
                   end
               end
        else
            for a in 1:cols
                Hx[rows, a] = (a % 4 == 0) || (a % 4 == 1) ? 0 : 1
            end
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
