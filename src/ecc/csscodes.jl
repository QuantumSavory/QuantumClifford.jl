#IN PROGRESS
struct CSS <: AbstractECC

end
#structure for CSS codes: requires parameters/components -> classical parity check matrix 4 ex, qubit n 4 ex, ...
#HAVE A GENERATOR OF HAMMING CODES - used for testing

#----------------Defining classical code ------
#classical_code_1: generator matrix
classical_gm_input = [0:0 0:0 0:0]

#classical_code_2: parity-check matrix
classical_pm_input = [1:0 1:1 0:1]

#-----------Building dual code ----------------
#REVISE THIS MATHEMATICAL PROCESS 
#nemo library (doc)- remember : def space matrix (float -> binary)
#create matrix, over space range x, convert into nemo matrix, add simple
#=
MatrixSpace(ResidueRing(ZZ,2), rᴬ, Δ)

function dual_code(H)
    null = nullspace(H)[2]
    @assert all(a*null .== 0)
    @assert size(a,1) + size(null,2) == size(a,2)
    transpose(null)
end
=#

G2_orthcolumnspace(c::Rep3) = nullspace(M, rtol=3) 

GD(c::Rep3) = [0:G2 G2_orthcolumnspace:0]

#-----------Building CSS code ----------------

parity_checks(c::CSS) = S""

code_n(c::CSS) = #variable input dependent?

parity_matrix(c::CSS) = stab_to_gf2(parity_checks(c::CSS))

syndrome_circuit(c::CSS) = #TODO

#Enconding circuit ----------------------------------

encoding_circuit(c::CSS) = []#TODO
#----------------------------------------------------------------

code_k(c::CSS) = #TODO

code_s(c::CSS) = #TODO

rate(c::CSS) = #TODO

distance(c::CSS) = #TODO

logx_ops(c::CSS) = #TODO

logz_ops(c::CSS) = #TODO

logy_ops(c::CSS) = #TODO
