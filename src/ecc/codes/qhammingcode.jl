struct QHamming
    r::Int
end

struct InvalidRError <: Exception end  # Define a custom exception type

function parity_checks(c::QHamming)
    if !(3 <= c.r < 15)
        throw(InvalidRError())  # Throw the exception if r is out of bounds
    end
    n = 2^c.r
    k = 2^c.r - c.r - 2
   
    # Initialize parity check matrix
    H = zeros(Int, k, n)
   
    # Generate parity checks for the quantum Hamming code
    for i in 1:k
        for j in 1:n
            if j & (1 << (i - 1)) != 0
                H[i, j] = 1
            end
        end
    end
   
    return H
end

code_n(c::QHamming) = 2^c.r
