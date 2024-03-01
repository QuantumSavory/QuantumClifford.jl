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
