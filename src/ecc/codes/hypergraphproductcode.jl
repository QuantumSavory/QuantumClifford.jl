#The [[n,k,d]] CSS code's constuction is based on two binary linear seeds codes, i element of {1,2} C sub i with parameters [n sub i, k sub i, d sub i] defined as the kernel of r sub i x n sub i check matrices H sub i of rank n sub i - k sub i. The hypergraph product yields two classical codes C sub x,z with parity-check matrices

#H sub x = (H sub 1 ⊗ I sub n2   I sub r1 ⊗ (H sub 2)^T) (1)
#H sub z =( I sub n1 ⊗ H sub 2    (H sub 1)^T ⊗ I sub r2) (2)
#where I sub m is the m-dimensonal identity matrix. These two codes then yied a hypergraph product code via the CSS construction

struct HypergraphProduct <: AbstractECC
  n_i::Vector{Int}      # Vector of n_i values
  k_i::Vector{Int}      # Vector of k_i values
  d_i::Vector{Int}      # Vector of d_i values
  r_i::Vector{Int}      # Vector of r_i values
  Hx::Matrix{Bool}     # H_x parity-check matrix
  Hz::Matrix{Bool}     # H_z parity-check matrix

  function HypergraphProduct(n_i, k_i, d_i, r_i)
    if length(n_i) != 2 || length(k_i) != 2 || length(d_i) != 2 || length(r_i) != 2
      error("Hypergraph product requires exactly two seed codes (length of n_i, k_i, d_i, and r_i must be 2)")
    end
    # Check for valid seed code parameters
    for i in 1:2
      if !(n_i[i] in [1, 2]) || !(k_i[i] <= n_i[i])
        error("Seed code parameters for code $i are invalid (n_i[$i]: $(n_i[i]), k_i[$i]: $(k_i[i]))")
      end
    end
    new(n_i, k_i, d_i, r_i)
  end

   function compute_Hx(n_i, k_i, d_i, r_i)
    H1 = IMatrix(n_i[1]) - IMatrix(k_i[1])
    H2 = IMatrix(n_i[2]) - IMatrix(k_i[2])
    kron(H1, IMatrix(n_i[2])) * kron(IMatrix(r_i[1]), H2')
  end

  function compute_Hz(n_i, k_i, d_i, r_i)
    H1 = IMatrix(n_i[1]) - IMatrix(k_i[1])
    H2 = IMatrix(n_i[2]) - IMatrix(k_i[2])
    kron(IMatrix(n_i[1]), H2) * kron(H1', IMatrix(r_i[2]))
  end

  function code_n(c::HypergraphProduct)
    n_x = size(compute_Hx(c.n_i, c.k_i, c.d_i, c.r_i), 2)
    n_z = size(compute_Hz(c.n_i, c.k_i, c.d_i, c.r_i), 2)
    return n_x * n_z
  end

  function parity_checks_xz(c::HypergraphProduct)
    hx = compute_Hx(c.n_i, c.k_i, c.d_i, c.r_i)
    hz = compute_Hz(c.n_i, c.k_i, c.d_i, c.r_i)
    return hx[:, 1:end], hz[:, 1:end]  # Remove first column (redundant)
  end

  function parity_checks_x(c::HypergraphProduct)
    hx, _ = parity_checks_xz(c)
    return hx
  end

  function parity_checks_z(c::HypergraphProduct)
    _, hz = parity_checks_xz(c)
    return hz
  end

  # Full parity checks based on CSS construction
  function parity_checks(c::HypergraphProduct)
    hx, hz = parity_checks_xz(c)

    # Construct extended parity-check matrices
    extended_Hx = vcat(hx, zeros(size(hz)))
    extended_Hz = vcat(zeros(size(hx)), hz)

    # Create and return stabilizer representation
    Stabilizer(fill(0x0, size(hx, 1) + size(hz, 1)), extended_Hx, extended_Hz)
  end
end

# Helper function for identity matrix
function IMatrix(n)
  return eye(n, n)
end



#println("Full parity checks:")
#println(parity_checks(code))
