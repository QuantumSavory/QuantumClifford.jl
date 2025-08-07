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

const THROW_REPR_NEEDS_COMMUTATIVE_ALGEBRA = 
"The group algebra must be commutative when using a single `repr` function,\
which is not the case here. Please specify separate `A_repr` and `B_repr` \
instead of a single `repr`. The default choice of\
`A_repr=right_repr_matrix, B_repr=left_repr_matrix` is frequently sufficient."

const THROW_MUST_BE_DIVISIBLE_BY_3 =
"The input parameters must be divisible by 3"

const THROW_REQUIRED_POSITIVE_ARG = 
"The input parameter should be positive"

const THROW_INVALID_XZ_COMPONENTS = 
"`xzcomponents` should be `:split` or `:together`"

const THROW_INVALID_CODE_DIMENSION =
"Dimension of the input code must be at least 2"    

const THROW_INVALID_X_METACHECKS =
"`X`-metachecks (`Mx`) require `D ≥ 4`"    

const THROW_INVALID_Z_METACHECKS = 
"`Z`-metachecks (`Mz`) require `D ≥ 3`"    


const THROW_INVALID_CUDA_ARG = 
"First argument to @run_cuda should be a function call"

const THROW_MODEL_MEMORY_LIMIT =
"Model exceeded memory limits"

const THROW_MODEL_TIME_LIMIT =
"Model exceeded time limit"

const THROW_INVALID_PARAMETERS_GOLAY =
"Invalid parameters: `n` must be either 24 or 23 to obtain a valid code."

const THROW_INVALID_PARAMETERS_HAMMING =
"Invalid parameters: `r` must be ≥ 2 to obtain a valid code."

const THROW_INVALID_PARAMETERS_REEDMULLER =
"Invalid parameters: r must be non-negative and r ≤ m in order to valid code."
