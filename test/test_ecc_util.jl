using JuMP
using GLPK

# Computing minimum distance for quantum LDPC codes is a NP-Hard problem,
# but we can solve mixed integer linear program (MILP) for small instance sizes.
function minimum_distance(hx,lx)
    n = size(hx, 2) # number of qubits
    m = size(hx, 1) # number of stabilizers
    # Maximum stabilizer weight
    whx = maximum(sum(hx[i, :] for i in 1:m))  # Maximum row sum of hx
    # Weight of the logical operator
    wlx = count(!iszero, lx)  # Count non-zero elements in lx
    # Number of slack variables needed to express orthogonality constraints modulo two
    num_anc_hx = ceil(Int, log2(whx))
    num_anc_logical = ceil(Int, log2(wlx))
    num_var = n + m * num_anc_hx + num_anc_logical
    model = Model(GLPK.Optimizer) # model
    set_silent(model)
    @variable(model, x[1:num_var], Bin) # binary variables
    @objective(model, Min, sum(x[i] for i in 1:n)) # objective function
    # Orthogonality to rows of hx constraints
    for row in 1:m
        weight = zeros(Int, num_var)
        supp = findall(hx[row, :] .!= 0)  # Non-zero indices in hx[row, :]
             for q in supp
                 weight[q] = 1
             end
             cnt = 1
             for q in 1:num_anc_hx
                 weight[n + (row - 1) * num_anc_hx + q] = -(1 << cnt)
             cnt += 1
         end
        @constraint(model, sum(weight[i] * x[i] for i in 1:num_var) == 0)
    end
    # Odd overlap with lx constraint
    supp = findall(lx .!= 0)  # Non-zero indices in lx
    weight = zeros(Int, num_var)
    for q in supp
        weight[q] = 1
    end
    cnt = 1
    for q in 1:num_anc_logical
        weight[n + m * num_anc_hx + q] = -(1 << cnt)
        cnt += 1
    end
    @constraint(model, sum(weight[i] * x[i] for i in 1:num_var) == 1)
    optimize!(model)
    opt_val = sum(value(x[i]) for i in 1:n)
    return Int(opt_val)
end
