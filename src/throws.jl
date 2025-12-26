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

const THROW_REINTERPRET_SIZE_MISMATCH = 
"Unable to perform the requested operation as the bit count of the principal \
axis of the underlying data is not divisible by that of the target type."
