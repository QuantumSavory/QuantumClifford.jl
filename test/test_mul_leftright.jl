@testitem "Mul leftright" begin
  @testset "Pauli string multiplication" begin
    using QuantumClifford: mul_left!, mul_right!
    test_sizes = [1,2,10,63,64,65,127,128,129] # Including sizes that would test off-by-one errors in the bit encoding.


    for n in test_sizes
      for _ in 1:20
        p1 = random_pauli(n)
        p2 = random_pauli(n)
        s = random_stabilizer(n)
        i = rand(1:n)
        @test p1*p2 == mul_left!(copy(p2), p1)
        @test p1*p2 == mul_right!(copy(p1), p2)
        @test mul_left!(copy(p2), p1) == (-1)^comm(p1,p2) * mul_right!(copy(p2), p1)
        @test mul_left!(copy(p2), s[i]) == mul_left!(copy(p2), s, i) == s[i]*p2
        @test mul_right!(copy(p2), s[i]) == mul_right!(copy(p2), s, i) == p2*s[i]
        @test mul_left!(copy(s), p2)[i] == p2*s[i]
        @test mul_right!(copy(s), p2)[i] == s[i]*p2
      end
    end
  end

  # test for #320
  @testset "verify SIMD implementation" begin
    for i in 1:10,
      n in 1:30,
      T in [UInt8, UInt16, UInt32, UInt64]
      a = rand(T, n)
      b = rand(T, n)
      c1,c2 = QuantumClifford.mul_ordered!(copy(a),copy(b))
      n1,n2 = QuantumClifford._mul_ordered_nonvec!(copy(a),copy(b))
      np = ((n1 ⊻ (n2<<1))&0x3)
      cp = ((c1 ⊻ (c2<<1))&0x3)
      @test np==cp
    end
  end

  @testset "OpenCL" begin
    using pocl_jll, OpenCL
    using KernelAbstractions: synchronize, get_backend
    using QuantumClifford: mul_left!, mul_right!, Tableau
    # Small sizes for encoding issues, large sizes for race conditions.
    test_sizes = [31, 32, 33, 63, 64, 65, 127, 128, 129,
                  64 * 1023, 64 * 1024, 64 * 1025,
                  64 * 2047, 64 * 2048, 64 * 2049]

    skip_flag = true
    if !skip_flag
      for n in test_sizes
        for _ in 1:5
          h_p1 = random_pauli(n)
          d_p1 = PauliOperator(CLArray(h_p1.phase), h_p1.nqubits, CLArray(h_p1.xz))
          h_p2 = random_pauli(n)
          d_p2 = PauliOperator(CLArray(h_p2.phase), h_p2.nqubits, CLArray(h_p2.xz))
          h_s = random_stabilizer(n)
          d_s = Stabilizer(Tableau(
            CLArray(h_s.tab.phases), h_s.tab.nqubits, CLArray(h_s.tab.xzs)
            ))
          i = rand(1:n)
          backend = get_backend(d_p1.phase)

          d_o = mul_left!(copy(d_p2), d_p1)
          h_o = mul_left!(copy(h_p2), h_p1)
          synchronize(backend)
          @test h_o.phase == Array(d_o.phase)
          @test h_o.xz == Array(d_o.xz)

          d_o = mul_right!(copy(d_p1), d_p2)
          h_o = mul_right!(copy(h_p1), h_p2)
          synchronize(backend)
          @test h_o.phase == Array(d_o.phase)
          @test h_o.xz == Array(d_o.xz)

          d_L = mul_left!(copy(d_p2), d_p1)
          d_R = mul_right!(copy(d_p2), d_p1)
          synchronize(backend)
          @test all(Array(d_L.phase) .== (-1)^comm(h_p1, h_p2) .* Array(d_R.phase))
          @test Array(d_L.xz) == Array(d_R.xz)

          d_o = mul_left!(copy(d_p2), d_s, i)
          h_o = mul_left!(copy(h_p2), h_s, i)
          synchronize(backend)
          @test h_o.phase == Array(d_o.phase)
          @test h_o.xz == Array(d_o.xz)

          d_o = mul_right!(copy(d_p2), d_s, i)
          h_o = mul_right!(copy(h_p2), h_s, i)
          synchronize(backend)
          @test h_o.phase == Array(d_o.phase)
          @test h_o.xz == Array(d_o.xz)

          d_o = mul_left!(copy(d_s), d_p2)
          h_o = mul_left!(copy(h_s), h_p2)
          synchronize(backend)
          @test h_o.phases[i] == Array(d_o.phases)[i]
          @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])

          d_o = mul_right!(copy(d_s), d_p2)
          h_o = mul_right!(copy(h_s), h_p2)
          synchronize(backend)
          @test h_o.phases[i] == Array(d_o.phases)[i]
          @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])

          # Potential race condition.
          d_o = mul_left!(copy(d_s), i, i)
          h_o = mul_left!(copy(h_s), i, i)
          synchronize(backend)
          @test h_o.phases[i] == Array(d_o.phases)[i]
          @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])
        end
      end
    end
  end

end

@testitem "Mul leftright -- GPU" tags = [:gpu] begin
  @testset "CUDA" begin
    using CUDA
    if length(CUDA.devices()) > 0

      using KernelAbstractions: synchronize, get_backend
      using QuantumClifford: mul_left!, mul_right!, Tableau
      # Small sizes for encoding issues, large sizes for race conditions.
      test_sizes = [31, 32, 33, 63, 64, 65, 127, 128, 129,
                    64 * 1023, 64 * 1024, 64 * 1025,
                    64 * 2047, 64 * 2048, 64 * 2049]

      skip_flag = true
      if !skip_flag
        for n in test_sizes
          for _ in 1:5
            h_p1 = random_pauli(n)
            d_p1 = PauliOperator(CuArray(h_p1.phase), h_p1.nqubits, CuArray(h_p1.xz))
            h_p2 = random_pauli(n)
            d_p2 = PauliOperator(CuArray(h_p2.phase), h_p2.nqubits, CuArray(h_p2.xz))
            h_s = random_stabilizer(n)
            d_s = Stabilizer(Tableau(
              CuArray(h_s.tab.phases), h_s.tab.nqubits, CuArray(h_s.tab.xzs)
              ))
            i = rand(1:n)
            backend = get_backend(d_p1.phase)

            d_o = mul_left!(copy(d_p2), d_p1)
            h_o = mul_left!(copy(h_p2), h_p1)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_o = mul_right!(copy(d_p1), d_p2)
            h_o = mul_right!(copy(h_p1), h_p2)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_L = mul_left!(copy(d_p2), d_p1)
            d_R = mul_right!(copy(d_p2), d_p1)
            synchronize(backend)
            @test all(Array(d_L.phase) .== (-1)^comm(h_p1, h_p2) .* Array(d_R.phase))
            @test Array(d_L.xz) == Array(d_R.xz)

            d_o = mul_left!(copy(d_p2), d_s, i)
            h_o = mul_left!(copy(h_p2), h_s, i)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_o = mul_right!(copy(d_p2), d_s, i)
            h_o = mul_right!(copy(h_p2), h_s, i)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_o = mul_left!(copy(d_s), d_p2)
            h_o = mul_left!(copy(h_s), h_p2)
            synchronize(backend)
            @test h_o.phases[i] == Array(d_o.phases)[i]
            @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])

            d_o = mul_right!(copy(d_s), d_p2)
            h_o = mul_right!(copy(h_s), h_p2)
            synchronize(backend)
            @test h_o.phases[i] == Array(d_o.phases)[i]
            @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])

            # Potential race condition.
            d_o = mul_left!(copy(d_s), i, i)
            h_o = mul_left!(copy(h_s), i, i)
            synchronize(backend)
            @test h_o.phases[i] == Array(d_o.phases)[i]
            @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])
          end
        end
      end
    end
  end

  @testset "ROCm" tags = [:gpu] begin
    using AMDGPU
    if length(AMDGPU.devices()) > 0

      using KernelAbstractions: synchronize, get_backend
      using QuantumClifford: mul_left!, mul_right!, Tableau
      # Small sizes for encoding issues, large sizes for race conditions.
      test_sizes = [31, 32, 33, 63, 64, 65, 127, 128, 129,
                    64 * 1023, 64 * 1024, 64 * 1025,
                    64 * 2047, 64 * 2048, 64 * 2049]

      skip_flag = true
      if !skip_flag
        for n in test_sizes
          for _ in 1:5
            h_p1 = random_pauli(n)
            d_p1 = PauliOperator(ROCArray(h_p1.phase), h_p1.nqubits, ROCArray(h_p1.xz))
            h_p2 = random_pauli(n)
            d_p2 = PauliOperator(ROCArray(h_p2.phase), h_p2.nqubits, ROCArray(h_p2.xz))
            h_s = random_stabilizer(n)
            d_s = Stabilizer(Tableau(
              ROCArray(h_s.tab.phases), h_s.tab.nqubits, ROCArray(h_s.tab.xzs)
              ))
            i = rand(1:n)
            backend = get_backend(d_p1.phase)

            d_o = mul_left!(copy(d_p2), d_p1)
            h_o = mul_left!(copy(h_p2), h_p1)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_o = mul_right!(copy(d_p1), d_p2)
            h_o = mul_right!(copy(h_p1), h_p2)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_L = mul_left!(copy(d_p2), d_p1)
            d_R = mul_right!(copy(d_p2), d_p1)
            synchronize(backend)
            @test all(Array(d_L.phase) .== (-1)^comm(h_p1, h_p2) .* Array(d_R.phase))
            @test Array(d_L.xz) == Array(d_R.xz)

            d_o = mul_left!(copy(d_p2), d_s, i)
            h_o = mul_left!(copy(h_p2), h_s, i)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_o = mul_right!(copy(d_p2), d_s, i)
            h_o = mul_right!(copy(h_p2), h_s, i)
            synchronize(backend)
            @test h_o.phase == Array(d_o.phase)
            @test h_o.xz == Array(d_o.xz)

            d_o = mul_left!(copy(d_s), d_p2)
            h_o = mul_left!(copy(h_s), h_p2)
            synchronize(backend)
            @test h_o.phases[i] == Array(d_o.phases)[i]
            @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])

            d_o = mul_right!(copy(d_s), d_p2)
            h_o = mul_right!(copy(h_s), h_p2)
            synchronize(backend)
            @test h_o.phases[i] == Array(d_o.phases)[i]
            @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])

            # Potential race condition.
            d_o = mul_left!(copy(d_s), i, i)
            h_o = mul_left!(copy(h_s), i, i)
            synchronize(backend)
            @test h_o.phases[i] == Array(d_o.phases)[i]
            @test (@view h_o.xzs[:, i]) == (@view Array(d_o.xzs)[:, i])
          end
        end
      end
    end
  end

end
