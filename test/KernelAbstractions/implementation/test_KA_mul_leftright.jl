# TODO: Revisit should the base package support more multiplication signatures.
# This must be done explicitly as they are not exported.
using QuantumClifford: mul_left!, mul_right!

@inline function test_KA_mul_leftright(AT, synchronize)
    cache = AllocCache()
    for mul! in (mul_left!, mul_right!)
    for n in test_sizes
    # Keep the memory usage sane.
    rows = min(n, max_rows)
    for r in one(round_count) : round_count
    block = block_sizes[r]
    batch = batch_sizes[r]
    @info "Round: Pauli size = $n, block size = $block, batch size = $batch."
    @cached cache begin

        # PauliOperator
        h_p1 = random_pauli(n)
        d_p1 = PauliOperator(AT(u32(h_p1.phase)), h_p1.nqubits, AT(h_p1.xz))
        h_p2 = random_pauli(n)
        d_p2 = PauliOperator(AT(u32(h_p2.phase)), h_p2.nqubits, AT(h_p2.xz))

        # Stabilizer
        h_s = Stabilizer(
            Tableau(
                rand(eltype(h_p1.phase), rows) .& 0x3,
                n,
                rand(eltype(h_p1.xz), (length(h_p1.xz), rows))
                )
            )
        d_s = Stabilizer(
            Tableau(
                AT(u32(h_s.tab.phases)),
                h_s.tab.nqubits,
                AT(h_s.tab.xzs)
                )
            )

        # Destabilizer
        h_d = Destabilizer(
            Tableau(
                rand(eltype(h_p1.phase), rows << 1) .& 0x3,
                rows,
                rand(
                    eltype(h_p1.xz),
                    cld(rows << 1, count_zeros(zero(eltype(h_p1.xz)))),
                    rows << 1
                    )
                )
            )
        d_d = Destabilizer(
            Tableau(
                AT(u32(h_d.tab.phases)),
                h_d.tab.nqubits,
                AT(h_d.tab.xzs)
                )
            )
        i = rand(one(rows) : rows)
        j = rand(one(rows) : rows)

        # Left/Right order.
        d_L = mul_left!(
            copy(d_p1), d_p2;
            block_size = Val(block), batch_size = Val(batch)
            )
        d_R = mul_right!(
            copy(d_p1), d_p2;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        # Either commutes or anti-commutes.
        @test begin
            all(
                (Array(d_L.phase) .- Array(d_R.phase)) .& 0x3
                .== (comm(h_p1, h_p2) << 1)
                )
            Array(d_L.xz) == Array(d_R.xz)
        end

        # Potential aliasing problem.
        h_o = mul!(copy(get_pauli(h_s, i)), h_s, i)
        d_o = mul!(
            copy(d_s), i, i;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            h_o.phase == Array(phases(d_o, i))
            h_o.xz == Array(xzs(d_o, i))
        end

        # PauliOperator - PauliOperator
        h_o = mul!(copy(h_p1), h_p2)
        d_o = mul!(
            copy(d_p1), d_p2;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            h_o.phase == Array(d_o.phase)
            h_o.xz == Array(d_o.xz)
        end

        for (h_v, d_v) in ((h_s.tab, d_s.tab), (h_s, d_s))

        # PauliOperator - Tableau/AbstractStabilizer[i]
        h_o = mul!(copy(h_p1), h_v, i)
        d_o = mul!(
            copy(d_p1), d_v, i;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            h_o.phase == Array(d_o.phase)
            h_o.xz == Array(d_o.xz)
        end

        # Marks the end for (h_v, d_v)
        end

        for (h_u, d_u) in ((h_s.tab, d_s.tab), (h_s, d_s))

        # Tableau/AbstractStabilizer - PauliOperator
        h_o = mul!(copy(h_u), h_p1)
        d_o = mul!(
            copy(d_u), d_p1;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            phases(h_o, i) == Array(phases(d_o, i))
            xzs(h_o, i) == Array(xzs(d_o, i))
        end

        # Tableau/AbstractStabilizer[i] - PauliOperator
        h_o = mul!(copy(get_pauli(h_u, i)), h_p1)
        d_o = mul!(
            copy(d_u), i, d_p1;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            h_o.phase == Array(phases(d_o, i))
            h_o.xz == Array(xzs(d_o, i))
        end

        # Tableau/AbstractStabilizer[i] - Self[j]
        h_o = mul!(copy(get_pauli(h_u, i)), h_u, j)
        d_o = mul!(
            copy(d_u), i, j;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            h_o.phase == Array(phases(d_o, i))
            h_o.xz == Array(xzs(d_o, i))
        end

        for (h_v, d_v) in ((h_s.tab, d_s.tab), (h_s, d_s))

        # Tableau/AbstractStabilizer - Tableau/AbstractStabilizer
        d_o = mul!(
            copy(d_u), d_v;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        # Rows become {+/-} x Identity since it is multiplied by itself.
        @test begin
            reduce(|, phases(d_o)) & 0x1 == 0
            reduce(|, xzs(d_o, i)) == 0
        end

        # Tableau/AbstractStabilizer - Tableau/AbstractStabilizer[i]
        h_o = mul!(copy(h_u), get_pauli(h_v, i))
        d_o = mul!(
            copy(d_u), d_v, i;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            phases(h_o, j) == Array(phases(d_o, j))
            xzs(h_o, j) == Array(xzs(d_o, j))
        end

        # Tableau/AbstractStabilizer[i] - Tableau/AbstractStabilizer[j]
        h_o = mul!(copy(get_pauli(h_u, i)), h_v, j)
        d_o = mul!(
            copy(d_u), i, d_v, j;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            h_o.phase == Array(phases(d_o, i))
            h_o.xz == Array(xzs(d_o, i))
        end

        # Marks the end for (h_v, d_v)
        end

        # Marks the end for (h_u, d_u)
        end

        # (Mixed)Destabilizer[i] - Self[j]
        h_o1 = mul!(
            copy(get_pauli(h_d, j)), get_pauli(h_d, i);
            phases = Val(false)
            )
        n_q = h_d.tab.nqubits
        h_o2 = mul!(copy(get_pauli(h_d, i + n_q)), get_pauli(h_d, j + n_q))
        d_o = mul!(
            d_d, i, j;
            block_size = Val(block), batch_size = Val(batch)
            )
        synchronize()
        @test begin
            h_o1.phase == Array(phases(d_o, j))
            h_o1.xz == Array(xzs(d_o, j))
            h_o2.phase == Array(phases(d_o, i + n_q))
            h_o2.xz == Array(xzs(d_o, i + n_q))
        end

    # Marks the end for @cached
    end
    # Marks the end for r in one(round_count) : round_count
    end
    # Marks the end for n in test_sizes
    end
    # Marks the end for mul!
    end
    unsafe_free!(cache)
end
