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
	parent_xzs = t.xzs isa LinearAlgebra.Adjoint ? parent(t.xzs) : t.xzs isa LinearAlgebra.Transpose ? parent(t.xzs) : t.xzs
	old = eltype(parent_xzs)
	chunk_axis = size(parent_xzs, 1)
	rows_axis = size(parent_xzs, 2)
	old_chunks_expected = 2 * QuantumClifford._nchunks(t.nqubits, old)
	if chunk_axis != old_chunks_expected
		total_bytes = length(parent_xzs) * sizeof(old)
		q, r = divrem(total_bytes, sizeof(U))
		throw(ArgumentError(THROW_REINTERPRET_SIZE_MISMATCH(:tableau, old, U, total_bytes, q, r)))
	end
	new_chunk_axis = 2 * QuantumClifford._nchunks(t.nqubits, U)
	if iszero(new_chunk_axis) || isodd(new_chunk_axis)
		total_bytes = length(parent_xzs) * sizeof(old)
		q, r = divrem(total_bytes, sizeof(U))
		throw(ArgumentError(THROW_REINTERPRET_SIZE_MISMATCH(:tableau, old, U, total_bytes, q, r)))
	end
	total_bytes = length(parent_xzs) * sizeof(old)
	q, r = divrem(total_bytes, sizeof(U))
	if r != 0 || q != length(parent_xzs) * sizeof(old) ÷ sizeof(U)
		throw(ArgumentError(THROW_REINTERPRET_SIZE_MISMATCH(:tableau, old, U, total_bytes, q, r)))
	end

	reint = reinterpret(U, vec(parent_xzs))
	new_mat = reshape(reint, new_chunk_axis, rows_axis)
	new_mat = copy(new_mat) # materialize to plain Array if needed

	if t.xzs isa LinearAlgebra.Adjoint
		return Tableau(t.phases, t.nqubits, new_mat')
	elseif t.xzs isa LinearAlgebra.Transpose
		return Tableau(t.phases, t.nqubits, transpose(new_mat))
	else
		return Tableau(t.phases, t.nqubits, new_mat)
	end
end

function Base.reinterpret(::Type{U}, s::Stabilizer) where {U <: Unsigned}
	try
		return Stabilizer(reinterpret(U, tab(s)))
	catch e
		if isa(e, ArgumentError) && occursin("Unable to reinterpret tableau storage", sprint(showerror, e))
			# Re-throw with stabilizer context but preserve the original error message
			err_msg = replace(sprint(showerror, e), "tableau" => "stabilizer")
			throw(ArgumentError(err_msg))
		end
		rethrow(e)
	end
end

function Base.reinterpret(::Type{U}, s::Destabilizer) where {U <: Unsigned}
	try
		return Destabilizer(reinterpret(U, tab(s)))
	catch e
		if isa(e, ArgumentError) && occursin("Unable to reinterpret tableau storage", sprint(showerror, e))
			# Re-throw with destabilizer context but preserve the original error message
			err_msg = replace(sprint(showerror, e), "tableau" => "destabilizer")
			throw(ArgumentError(err_msg))
		end
		rethrow(e)
	end
end

function Base.reinterpret(::Type{U}, ms::MixedStabilizer) where {U <: Unsigned}
	try
		return MixedStabilizer(reinterpret(U, tab(ms)), rank(ms))
	catch e
		if isa(e, ArgumentError) && occursin("Unable to reinterpret tableau storage", sprint(showerror, e))
			# Re-throw with mixedstabilizer context but preserve the original error message
			err_msg = replace(sprint(showerror, e), "tableau" => "mixedstabilizer")
			throw(ArgumentError(err_msg))
		end
		rethrow(e)
	end
end

function Base.reinterpret(::Type{U}, ms::MixedDestabilizer) where {U <: Unsigned}
	try
		return MixedDestabilizer(reinterpret(U, tab(ms)), rank(ms))
	catch e
		if isa(e, ArgumentError) && occursin("Unable to reinterpret tableau storage", sprint(showerror, e))
			# Re-throw with mixeddestabilizer context but preserve the original error message
			err_msg = replace(sprint(showerror, e), "tableau" => "mixeddestabilizer")
			throw(ArgumentError(err_msg))
		end
		rethrow(e)
	end
end

function Base.reinterpret(::Type{U}, f::PauliFrame) where {U <: Unsigned}
	try
		reinterpreted_frame = reinterpret(U, f.frame)
		new_tab = Tableau(copy(reinterpreted_frame.tab.phases), reinterpreted_frame.tab.nqubits, collect(reinterpreted_frame.tab.xzs))
		return PauliFrame(typeof(f.frame)(new_tab), f.measurements)
	catch e
		if isa(e, ArgumentError) && (occursin("Unable to reinterpret tableau storage", sprint(showerror, e)) || occursin("Unable to reinterpret stabilizer storage", sprint(showerror, e)))
			# Re-throw with pauliframe context but preserve the original error message
			err_msg = replace(sprint(showerror, e), "tableau" => "pauliframe", "stabilizer" => "pauliframe")
			throw(ArgumentError(err_msg))
		end
		rethrow(e)
	end
end
