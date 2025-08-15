
#=============================================================================#
# There is vanishing overhead for supporting both of these.
@enum MultiplicationOrder::UInt8 begin
    multiplication_order_left
    multiplication_order_right
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
    sort_order_pauli_bit_prefer_x
    sort_order_pauli_bit_prefer_z
    sort_order_qubit_number_prefer_x
    sort_order_qubit_number_prefer_z
end

# Determines whether contraction proceeds from high to low or low to high.
@enum ScanSide::UInt8 begin
    scan_side_lesser
    scan_side_greater
end

# Provides enhanced clarity over plain numerical values.
# CAUTION: The values are NOT arbitrary but utilised for proper indexing.
@enum TrackerContent::UInt8 begin
    tracker_content_index = 0x1
    tracker_content_bit_shift = 0x2
    tracker_content_bit_type = 0x3
    tracker_content_swap_from = 0x4
    tracker_content_swap_to = 0x5
end

# Provides enhanced clarity over plain numerical values.
# CAUTION: The values are NOT arbitrary but utilised for proper ordering.
@enum PauliBit::UInt8 begin
    pauli_bit_x = 0x0
    pauli_bit_z = 0x1
    pauli_bit_invalid = 0x2
end

# Provides enhanced clarity over plain numerical values.
# CAUTION: The values are NOT arbitrary but utilised for proper masking.
@enum BitFieldFlags::DeviceUnsigned begin
    bit_field_flag_leader = 0x1
    bit_field_flag_scan = 0x2
    bit_field_flag_multiply = 0x4
end
#=============================================================================#
