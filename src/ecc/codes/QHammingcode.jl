#The Gottesman code has the parameters [[2^r, 2^r - r - 2, 3]]
#So, the parity check matrix will be of the form:
#[[I_r, K, -I_r]
#[K^T, I_(2^r - r - 2), 0]
#[0, -K, I_r]]

struct QHamming
  r::Int
end

function parity_checks(code::QHamming)
  r = code.r
  n = 2^r

  # Construct the K matrix (refer to Gottesman-Chuang codes for details)
  K = zeros(Int, (n - r - 2), r)
  for i in 1:n-r-2
    j = i + 1
    while j <= n
      K[i, bit2int(Integer.tobits(j - 1)[end])] = 1
      j += 1
    end
  end

  # Parity check matrix
  H = [[diagm(ones(Int, r)), K, -diagm(ones(Int, r))]
       [K', diagm(ones(Int, n - r - 2)), zeros(Int, n - r - 2)]
       [zeros(Int, n - r - 2), -K, diagm(ones(Int, r))]]

  return CSS(H)
end

#how to use
# Create two different instances of the quantum Hamming code
#code1 = QHamming(5)  # r = 5
#code2 = QHamming(7)  # r = 7

# Get the parity check matrices for each instance
#H1 = parity_checks(code1)
#H2 = parity_checks(code2)

