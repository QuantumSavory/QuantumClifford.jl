
#=============================================================================#
# Dictates the direction of the multiplication operation(s).
@enum MultiplicationOrder::UInt8 begin
    multiplication_order_left
    multiplication_order_right
end

# Specifies whether ordering should prioritise X or Z pauli operators.
@enum PauliPreferance::UInt8 begin
    pauli_preferance_x
    pauli_preferance_z
end

# Determines whether the first grid dimension matches the rows or the qubits.
@enum PrimaryAxis::UInt8 begin
    primary_axis_rows
    primary_axis_qubits
end

#==============================================================================
UTILISED INTERNALLY
==============================================================================#

# Pauli bit: Highest X(Z) < Lowest Z(X). Qubit Number: 1X(Z) < 1Z(X) < 2X(Z)...
@enum SortOrder::UInt8 begin
    sort_order_pauli_bit
    sort_order_qubit_number
end

# Determines whether contraction proceeds from high to low or low to high.
@enum ScanSide::UInt8 begin
    scan_side_lesser
    scan_side_greater
end

# Provides enhanced clarity over plain numerical values.
# CAUTION: The values are NOT arbitrary but utilised for proper ordering.
@enum PauliBit::UInt8 begin
    pauli_bit_primary = 0
    pauli_bit_secondary = 1
    pauli_bit_invalid = 2
end

# Provides enhanced clarity over plain numerical values.
# CAUTION: The values are NOT arbitrary but utilised for proper masking.
@enum BitFieldFlags::DeviceUnsigned begin
    bit_field_flag_leader = 1
    bit_field_flag_scan = 2
    bit_field_flag_multiply = 4
end

# Provides enhanced clarity over plain numerical values.
# CAUTION: The values are NOT arbitrary but utilised for proper indexing.
@enum TrackerContent::Csize_t begin
    tracker_content_index = 1
    tracker_content_bit_shift = 2
    tracker_content_bit_type = 3
    tracker_content_swap_from = 4
    tracker_content_swap_to = 5
end
#=============================================================================#
