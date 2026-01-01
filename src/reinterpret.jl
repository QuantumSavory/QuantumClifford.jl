function Base.reinterpret(::Type{U}, p::PauliOperator) where {U <: Unsigned}
	old = eltype(p.xz)
	total_bytes = length(p.xz) * sizeof(old)
	q, r = divrem(total_bytes, sizeof(U))
	if r != 0 || isodd(q)
		throw(ArgumentError(THROW_REINTERPRET_SIZE_MISMATCH(:pauli, old, U, total_bytes, q, r)))
	end
	new_xz = reinterpret(U, p.xz)
	return PauliOperator(p.phase, p.nqubits, new_xz)
end

function Base.reinterpret(::Type{U}, t::Tableau) where {U <: Unsigned}
	rows = size(t.xzs, 1)
	old = eltype(t.xzs)
	total_bytes = rows * sizeof(old)
	q, r = divrem(total_bytes, sizeof(U))
	if r != 0 || isodd(q)
		throw(ArgumentError(THROW_REINTERPRET_SIZE_MISMATCH(:tableau, old, U, total_bytes, q, r)))
	end
	new_xzs = reinterpret(U, t.xzs)
	return Tableau(t.phases, t.nqubits, new_xzs)
end

Base.reinterpret(::Type{U}, s::Union{Stabilizer, Destabilizer}) where {U <: Unsigned} =
	typeof(s)(reinterpret(U, tab(s)))

Base.reinterpret(::Type{U}, ms::Union{MixedStabilizer, MixedDestabilizer}) where {U <: Unsigned} =
	typeof(ms)(reinterpret(U, tab(ms)), rank(ms))

function Base.reinterpret(::Type{U}, f::PauliFrame) where {U <: Unsigned}
	return PauliFrame(reinterpret(U, f.frame), f.measurements)
end
