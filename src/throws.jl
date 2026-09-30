const NoZeroQubit = ArgumentError("Qubit indices have to be larger than zero, but you are attempting to create a gate acting on a qubit with a non-positive index. Ensure indexing always starts from 1.")

const THROW_BOUNDS =
"Unable to perform the requested operation due to encountering a mismatch \
between the provided index(indices) and the pertinent size(s) of the other \
argument(s)."

const THROW_SIZE =
"Unable to perform the requested operation due to encountering a mismatch \
between the pertinent size(s) of the provided arguments."

const THROW_NQUBITS =
"Unable to perform the requested operation due to encountering a mismatch \
between the number of qubits in the provided arguments."

function THROW_REINTERPRET_SIZE_MISMATCH(kind::Symbol, oldT::Type, newT::Type,
		total_bytes::Integer, new_words::Integer, remainder::Integer)
	size_old = sizeof(oldT)
	size_new = sizeof(newT)
	prefix = "Unable to reinterpret $kind storage"

	align_hint = begin
		if kind === :pauli
			"PauliOperator X/Z halves"
		elseif kind === :tableau
			"Tableau X/Z row pairs"
		elseif kind === :stabilizer
			"Stabilizer X/Z row pairs"
		elseif kind === :destabilizer
			"Destabilizer X/Z row pairs"
		elseif kind === :mixedstabilizer
			"MixedStabilizer X/Z row pairs"
		elseif kind === :mixeddestabilizer
			"MixedDestabilizer X/Z row pairs"
		elseif kind === :pauliframe
			"PauliFrame frame rows"
		else
			"tableau rows"
		end
	end

	if remainder != 0
		return "$prefix: $total_bytes bytes from elements of $oldT are not divisible by sizeof($newT)=$size_new."
	elseif isodd(new_words)
		return "$prefix: reinterpret would produce $new_words words of $newT, but an even count is required to align $align_hint."
	else
		return "$prefix: incompatible element sizes (sizeof($oldT)=$size_old, sizeof($newT)=$size_new)."
	end
end
