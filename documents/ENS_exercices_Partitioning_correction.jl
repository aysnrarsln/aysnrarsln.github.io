using JuMP
using Gurobi


prod_cost = 50
inventory_cost = 65
lot_cost = [1500,1875]
lot_size = [25, 25]
L = 1:2 # lot index
realizations = []
partitions = [[5, 95]] # initial partitions


function solve_partition(partitions)
    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "OutputFlag", 0)
    P = length(partitions)
    @variable(M,θ)
    @variable(M,x>=0)
    @variable(M,y[1:P,L], Bin)
    @constraint(M,[p in 1:P, i in 1:2], x+sum(lot_size[j]*y[p,j] for j in L) >= partitions[p][i])
    @constraint(M, [p in 1:P, i in 1:2], θ >= prod_cost*x + inventory_cost*(x+sum(lot_size[j]*y[p,j] for j in L)-partitions[p][i]) + sum(lot_cost[j]*y[p,j] for j in L))
    @objective(M, Min, θ)
    print(M)
    optimize!(M)
    return value(θ),value(x),value.(y)  
end 

function separate_inventory(LB, UB, xval, yval)
    sep_inv = Model(Gurobi.Optimizer)
    set_optimizer_attribute(sep_inv,"OutputFlag", 0)
    @variable(sep_inv, LB <= ξ <= UB)
    @objective(sep_inv, Min, xval+sum(lot_size[i]*yval[i] for i in L) - ξ)
    optimize!(sep_inv)
    return value(ξ)
end


function separate_objective(LB, UB, θval, xval, yval)
    sep_obj = Model(Gurobi.Optimizer)
    set_optimizer_attribute(sep_obj,"OutputFlag", 0)
    @variable(sep_obj, LB <= ξ <= UB)
    @objective(sep_obj, Min, θval - prod_cost*xval -  inventory_cost*(xval+sum(lot_size[i]*yval[i] for i in L)-ξ) - sum(lot_cost[i]*yval[i] for i in L))
    optimize!(sep_obj)
    return value(ξ)
end


function calculate_dual_bound(realizations)
    Dual = Model(Gurobi.Optimizer)
    S = length(realizations)
    set_optimizer_attribute(Dual,"OutputFlag", 0)
    @variable(Dual,θ)
    @variable(Dual,x>=0)
    @variable(Dual,y[1:S,L], Bin)
    @constraint(Dual,[s in 1:S], x+sum(lot_size[i]*y[s,i] for i in L) >= realizations[s])
    @constraint(Dual,[s in 1:S], θ >= prod_cost*x + inventory_cost*(x+sum(lot_size[i]*y[s,i] for i in L)-realizations[s])+sum(lot_cost[i]*y[s,i] for i in L))
    @objective(Dual, Min, θ)
    print(Dual)
    optimize!(Dual)
    return objective_value(Dual)
end

function create_partitions(realizations)
    sort!(realizations)
    
    partitions = []
    for i in 1:(length(realizations)-1)
        push!(partitions, [realizations[i], 0.5*(realizations[i]+realizations[i+1])])
        push!(partitions, [0.5*(realizations[i]+realizations[i+1]),realizations[i+1]])
    end

    return partitions
end


function apply_uncertainty_partition(partitions)
    θval, xval, yval, dual_bound = Inf, nothing, nothing, -Inf
    for i in 1:3
        θval, xval, yval = solve_partition(partitions)
        println("Partitioned solution:")
        println(θval)
        println(xval)
        println(yval)
        println("Primal bound:", θval)

        for p in 1:length(partitions)
            println("Critical scenarios")
            ξ = separate_inventory(partitions[p][1], partitions[p][2], xval, yval[p, :])
            push!(realizations, ξ)
            ξ = separate_objective(partitions[p][1], partitions[p][2], θval, xval, yval[p, :])
            push!(realizations, ξ)
            println(realizations)
        end

        unique_realizations = unique(realizations)
        dual_bound = calculate_dual_bound(unique_realizations)
        println("Dual bound:", dual_bound)

        partitions = create_partitions(unique_realizations)
        println("New partitions:", partitions)
    end
    return θval, xval, yval, dual_bound
end

θval, xval, yval, dual_bound = apply_uncertainty_partition(partitions)
println("Best primal bound obtained:", θval)
println("Value of x:", xval)  
println("Best dual bound obtained:", dual_bound)
println("Solution x=$xval has an optimality gap of ", (θval-dual_bound)/θval*100, "%")







