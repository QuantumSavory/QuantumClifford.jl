
#=============================================================================#
# CAUTION: Requires unsafe_indices = true if active_threads_count < block_size
# TODO: Overhaul once __ctx__ is no longer necessary for runtime queries.
# TODO: Revisit once warp level primitives are supported.
@inline @generated function shared_memory_reduce!(
    f!::Function, index::Integer, ::Val{block_size}, arguments...
    )::Nothing where {block_size}

    current = block_size
    if current > zero(current)
        stages = :(;;)
    end
    while current > one(current)
        stages = :(
            $stages;
            # The call to KA.@synchronize is ALL or NOTHING.
            KA.@synchronize();
            )
        if iseven(current)
            current >>= one(current)
            stages = :(
                $stages;
                if index <= $current;
                    f!(index, $current, arguments...);
                end;
                )
        else
            current = (current >> one(current)) + one(current)
            stages = :(
                $stages;
                # The strict inequality is intentional.
                if index < $current;
                    f!(index, $current, arguments...);
                end;
                )
        end
    end
    return :(
        $stages;
        return nothing;
        )

end

#==============================================================================
REDUCTION CATALOGUE
==============================================================================#

@inline @generated function reduce_sum!(
    index::Integer, stride::Integer, arguments::AbstractArray...
    )::Nothing

    if length(arguments) > 0x0
        reduction = :(
            arguments[0x1][index] += arguments[0x1][index + stride];
            )
    end
    for n in 0x2 : length(arguments)
        reduction = :(
            $reduction;
            arguments[$n][index] += arguments[$n][index + stride];
            )
    end
    return :(
        @inbounds begin;
            $reduction;
        end;
        return nothing;
        )

end

@inline @generated function reduce_lexicographic_min!(
    index::Integer, stride::Integer, arguments::AbstractArray...
    )::Nothing

    if length(arguments) > 0x0
        clause = :(arguments[0x1][index + stride] < arguments[0x1][index])
        body = :(
            arguments[0x1][index] = arguments[0x1][index + stride];
            )
    end
    for n in 0x2 : length(arguments)
        subclause = :(arguments[0x1][index + stride] == arguments[0x1][index])
        for m in 0x2 : (n - 0x1)
            subclause = :(
                $subclause &&
                    arguments[$m][index + stride] == arguments[$m][index]
                )
        end
        subclause = :(
            $subclause && arguments[$n][index + stride] < arguments[$n][index]
            )
        clause = :($clause || $subclause)
        body = :(
            $body;
            arguments[$n][index] = arguments[$n][index + stride];
            )
    end
    return :(
        @inbounds begin;
            if $clause;
                $body;
            end;
        end;
        return nothing;
        )

end
#=============================================================================#
