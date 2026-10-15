
#=============================================================================#
# CAUTION: Requires an integral ratio between the sizes of the data types.
@inline function change_type(
    ::Type{T}, source::AbstractArray{S}; interpret::Bool = false
    ) where {T <: Unsigned, S <: Unsigned}

    if T == S
        output = copy(source)
    elseif !interpret || size(source) == ()
        output = map(x -> T(x), source)
    elseif sizeof(S) >= sizeof(T)
        output = copy(reinterpret(T, source))
    else
        len, dims... = size(source)
        new_len = cld(len, div(sizeof(T), sizeof(S)))
        output = similar(source, T, new_len, dims...)
        @inbounds fill!(
            (@view output[new_len, Base.OneTo.(dims)...]),
            zero(T)
            )
        temp = reinterpret(S, output)
        @inbounds (@view temp[Base.OneTo(len), Base.OneTo.(dims)...]) .= source
    end
    return output

end
#=============================================================================#
