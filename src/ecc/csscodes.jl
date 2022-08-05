using Nemo 
using LinearAlgebra
#using Statistics

#structure for CSS codes
struct CSS <: AbstractECC end

#=struct CSS <: AbstractECC

    function css_n end
    function classical_code_G_matrix end
    function classical_code_H_matrix end
    function classic_parity_checks end
    function dual_code_prt1 end
    function dual_code_prt2 end

end=#

#=Testing code for CSS code contruction--------------------------------------------------
#[7,4] Hamming code -> Steane's 7

#classical_code_H_matrix
classical_code_H_matrix = [0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]
#classical_code_G_matrix - dual code
classical_code_G_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] #tst
#classical_code_G_matrix(c) = gf2_H_to_G(classical_code_H_matrix) 
---------------------------------------------------------------------------------------=#

include("hamming_code_generator.jl") #hamming code generator
import .Hamming
using .Hamming
using .Hamming: construct_A, construct_H_from_A, construct_G_from_A

A = construct_A()
classical_code_H_matrix = construct_H_from_A(A)
classical_code_G_matrix = construct_G_from_A(A)

#----------------------------------------------------------------------------------------------

#----------------CSS code generation ----------------------------------------------------------
#----------Dual code -------------------
#Size
size_row_H = size(classical_code_H_matrix, 1)
size_column_H = size(classical_code_H_matrix, 2)

size_row_G = size(classical_code_G_matrix, 1)
size_column_G = size(classical_code_G_matrix, 2)

#Dual code build
X_zeros = zeros(Int8, size_row_H, size_column_H)
Z_zeros = zeros(Int8, size_row_G, size_column_G)

#Final X & Z matrix
X_matrix = X_zeros
Z_matrix = classical_code_G_matrix
#TODO: overlap X & Z -> Y !!!
hcat(Z_matrix,Z_zeros)
hcat(X_matrix,classical_code_H_matrix)

#-----------Building CSS code ----------------

parity_checks(c::CSS) = Stabilizer(Z_matrix,X_matrix) #READ MANUAL 

code_n(c::CSS) = size(X_matrix, 1) #variable input dependant

parity_matrix(c::CSS) = stab_to_gf2(parity_checks) 

#Encoding circuit ----------------------------------

#encoding_circuit(c::CSS) = [] #TODO -> START SYNDROME CIRCUIT
#----------------------------------------------------------------


code_s(c::CSS) = length(parity_checks(c))

code_k(c::CSS) = code_n(c) - code_s(c)

rate(c::CSS) = code_k(c)/code_s(c)


logx_ops(c::CSS) = P"XXXXXXXXX"

logz_ops(cc::CSS) = P"ZZZZZZZZZ"

logy_ops(cc::CSS) = P"YYYYYYYYY" 

#naive_syndrome(c::CSS) #Syndrome circuit
