struct QHamming <:AbstractECC
    r::Int
end

function parity_checks(c::QHamming)
    r = c.r
    n = 2^r
    k = 2^r - r - 2  # Calculate number of checks

    # Initialize parity check matrix
    H = zeros(Bool, k, n)

    # Generate parity checks for the quantum Hamming code
    for i in 1:k
        for j in 1:n
            if j & (1 << (i - 1)) != 0
                H[i, j] = true
            end
        end
    end

   
    # Extract Hx and Hz from H
    Hx = H[:, 1:end-r]
    Hz = H[:, end-r+1:end]

    # Construct the Stabilizer object using the built-in constructor
    Stabilizer(fill(0x0, 2k), Hx, Hz)
end


