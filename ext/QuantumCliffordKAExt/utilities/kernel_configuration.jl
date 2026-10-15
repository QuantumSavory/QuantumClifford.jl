
#=============================================================================#
# Translates index set dimensions into the grid-block model.
@inline function tessellate(
    space::NTuple{N, Integer}, tile::NTuple{N, Integer}
    ) where {N}

    return tile .* cld.(space, tile)

end

# Commonly set unsafe_indices = true, hence replaces KA.@index(Global, NTuple).
@inline function global_index(
    block_index::NTuple{N, T},
    block_dim::NTuple{N, Integer},
    thread_index::NTuple{N, Integer}
    ) where {N, T <: Integer}

    return (block_index .- one(T)) .* block_dim .+ thread_index

end
#=============================================================================#
