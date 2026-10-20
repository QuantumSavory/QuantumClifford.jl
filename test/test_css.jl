@testitem "Local Clifford conversion to CSS" begin
    using Random
    using QuantumClifford
    using QuantumClifford: mul_left!, rowswap!

    rng = Xoshiro(0xc55)

    # Compare the signed groups, not a particular choice or ordering of generators.
    function canonical_basis(s)
        basis, _, r = canonicalize!(copy(s); ranks=true)
        return basis[1:r]
    end

    is_css_basis(s) = all(p -> !any(xbit(p)) || !any(zbit(p)), s)

    function apply_circuit!(s, circuit)
        for gate in circuit
            apply!(s, gate)
        end
        return s
    end

    function scramble_rows!(s)
        length(s) < 2 && return s
        for _ in 1:4length(s)
            i, j = randperm(rng, length(s))[1:2]
            rowswap!(s, i, j)
            mul_left!(s, i, j)
        end
        return s
    end

    function check_conversion(s)
        original = copy(s)
        result = lc_to_css(s)
        @test s == original
        @test result !== nothing
        isnothing(result) && return nothing
        css, circuit = result
        @test css isa Stabilizer
        @test circuit isa Vector{SingleQubitOperator}
        @test nqubits(css) == nqubits(s)
        @test length(css) == length(canonical_basis(s))
        @test is_css_basis(css)
        @test canonical_basis(apply_circuit!(copy(s), circuit)) == canonical_basis(css)
        @test canonical_basis(apply_circuit!(copy(css), inv.(reverse(circuit)))) == canonical_basis(s)
        return css, circuit
    end

    bell = S"XX
             ZZ"
    repetition = S"ZZI
                   IZZ"
    steane = S"XXXXIII
               XXIIXXI
               XIXIXIX
               ZZZZIII
               ZZIIZZI
               ZIZIZIZ"
    # Two disconnected factors, an unused qubit, and nontrivial generator signs.
    disconnected = S"-XXIIII
                      ZZIIII
                     -IIZZII
                      IIIZZI"

    @testset "CSS generators remain CSS after changing the generator basis" begin
        for original in (bell, repetition, steane, disconnected), _ in 1:8
            scrambled = scramble_rows!(copy(original))
            result = check_conversion(scrambled)
            isnothing(result) && continue
            css, circuit = result
            @test isempty(circuit)
            @test canonical_basis(css) == canonical_basis(original)
        end
        # Multiplying XX and ZZ introduces a minus sign: XX * ZZ = -YY.
        css, circuit = lc_to_css(S"XX
                                 -YY")
        @test isempty(circuit)
        @test canonical_basis(css) == canonical_basis(bell)
    end

    @testset "Recover locally rotated CSS codes" begin
        originals = [bell, repetition, steane, disconnected]
        # Disjoint X/Z checks commute; CNOTs preserve CSS form while spreading
        # their support. Include a register crossing the packed-word boundary.
        for n in (1, 2, 4, 8, 17, 65)
            r = rand(rng, 1:min(n, 6))
            css = zero(Stabilizer, r, n)
            for i in 1:r
                css[i, i] = (i <= r÷2, i > r÷2)
            end
            phases(css) .= rand(rng, UInt8[0, 2], r)
            if n > 1
                for _ in 1:4n
                    control, target = randperm(rng, n)[1:2]
                    apply!(css, sCNOT(control, target))
                end
            end
            @test is_css_basis(css)
            push!(originals, css)
        end

        for original in originals, _ in 1:8
            scrambled = scramble_rows!(copy(original))
            local_cliffords = [random_clifford1(rng, q) for q in 1:nqubits(original)]
            apply_circuit!(scrambled, local_cliffords)
            result = check_conversion(scrambled)
            isnothing(result) && continue
            css, circuit = result

            # CSS representatives are not unique. Undoing both known changes
            # must recover exactly the original signed CSS stabilizer group.
            recovered = apply_circuit!(copy(css), inv.(reverse(circuit)))
            apply_circuit!(recovered, inv.(reverse(local_cliffords)))
            @test canonical_basis(recovered) == canonical_basis(original)
        end
    end

    @testset "Second affine branch" begin
        # Fixing the first local b coordinate to zero is inconsistent here;
        # the alternative c-coordinate branch has a CSS conversion.
        check_conversion(S"YXXI
                           ZYIY")
    end

    @testset "Redundant generators and empty stabilizers" begin
        redundant = S"XX
                      ZZ
                     -YY
                      II
                      XX"
        for s in (redundant, zero(Stabilizer, 0, 0), zero(Stabilizer, 0, 5),
                  zero(Stabilizer, 3, 0), zero(Stabilizer, 3, 5))
            result = check_conversion(s)
            isnothing(result) && continue
            css, circuit = result
            @test isempty(circuit)
            @test canonical_basis(css) == canonical_basis(s)
        end
        rotated = apply_circuit!(copy(redundant), [sPhase(1), sHadamard(2)])
        check_conversion(rotated)
    end

    @testset "Reject non-CSS and invalid stabilizers" begin
        five_qubit = S"XZZXI
                       IXZZX
                       XIXZZ
                       ZXIXZ"
        original = copy(five_qubit)
        @test lc_to_css(five_qubit) === nothing
        @test five_qubit == original
        for invalid in (S"iX", S"X\nZ", S"-I", S"X\n-X", S"XX\nZZ\nYY")
            original = copy(invalid)
            @test_throws ArgumentError lc_to_css(invalid)
            @test invalid == original
        end
    end

    @testset "Agreement with exhaustive local Clifford search" begin
        # Only the six permutations of the Pauli axes matter for CSS existence.
        # This small oracle checks all local frames using existing row reduction,
        # independently of the linear-algebra criterion used by lc_to_css.
        function brute_force_lc_css(s)
            for choices in Iterators.product(ntuple(_ -> 1:6, nqubits(s))...)
                candidate = copy(s)
                for (q, choice) in enumerate(choices)
                    apply!(candidate, enumerate_single_qubit_gates(choice; qubit=q))
                end
                is_css_basis(canonical_basis(candidate)) && return true
            end
            return false
        end

        for n in 1:5, r in unique([1, max(1, n-1), n]), _ in 1:2
            s = random_stabilizer(rng, r, n)
            expected = brute_force_lc_css(s)
            result = lc_to_css(s)
            @test !isnothing(result) == expected
            expected && check_conversion(s)
        end
    end
end
