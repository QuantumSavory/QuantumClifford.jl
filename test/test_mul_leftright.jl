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
end
