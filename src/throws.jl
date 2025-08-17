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

const THROW_INVALID_XZ_COMPONENTS = 
"`xzcomponents` should be `:split` or `:together`"


const THROW_INVALID_PARAMETERS_REEDMULLER =
"Invalid parameters: r must be non-negative and r ≤ m in order to valid code."

const THROW_DELFOSSE_MIN_BLOCKS =
"The number of blocks must be at least 2 to construct a valid code."

const THROW_INVALID_ACTION_QUBITS =
"The tableau and the provided operator need to act on the same number of qubits.
Consider specifying an array of indices as a third argument to the `apply!` function to avoid this error."

const THROW_NO_ZERO_QUBIT = 
"Qubit indices have to be larger than zero, 
but you are attempting to create a gate acting on a qubit with a non-positive index. Ensure indexing always starts from 1."


