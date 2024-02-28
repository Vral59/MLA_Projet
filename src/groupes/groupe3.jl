# Guillaume Casanova et Raphael Taisant
# P-centres problems


function solve_pc1(n::Int, m::Int, cost_connection::Matrix{Int}, nb_centres::Int; relaxation::Bool=false, verbose::Int=0)
    """Solve the p-centers problem with formulation PC1"""
    # n number of clients
    # m number of sites

    model = Model(CPLEX.Optimizer)
    if verbose == 0
        set_silent(model)
    end

    @variable(model, x[1:n, 1:m] >= 0, Bin)
    @variable(model, y[1:m], Bin)
    @variable(model, r>=0)

    @constraint(model, [i in 1:n, j in 1:m], x[i,j] <= y[j])
    @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:m)  == 1)

    # p centres
    @constraint(model, sum(y[j] for j in 1:m) == nb_centres)

    # bound on r in PC1
    @constraint(model, [i in 1:n, j in 1:m], cost_connection[i,j] * x[i,j] <= r)

    @objective(model, Min, r)

    if relaxation
        relax_integrality(model)    
    end

    optimize!(model)

    return objective_value(model), solve_time(model), node_count(model), value.(x), value.(y)
end

function solve_pc(n::Int, m::Int, cost_connection::Matrix{Int}, nb_centres::Int; relaxation::Bool=false, verbose::Int=0)
    """Solve the p-centers problem with formulation PC1"""
    # n number of clients
    # m number of sites

    model = Model(CPLEX.Optimizer)
    if verbose == 0
        set_silent(model)
    end

    @variable(model, x[1:n, 1:m] >= 0, Bin)
    @variable(model, y[1:m], Bin)
    @variable(model, r>=0)

    @constraint(model, [i in 1:n, j in 1:m], x[i,j] <= y[j])
    @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:m)  == 1)

    # p centres
    @constraint(model, sum(y[j] for j in 1:m) == nb_centres)

    # bound on r in PC1
    @constraint(model, [i in 1:n], sum(cost_connection[i,j] * x[i,j] for j in 1:m)<= r)

    @objective(model, Min, r)

    if relaxation
        relax_integrality(model)    
    end

    optimize!(model)

    return objective_value(model), solve_time(model), node_count(model), value.(x), value.(y)
end

function lower_bound_0(n::Int, m::Int, cost_connection::Matrix{Int})
    """Compute a lower bound"""
    lb=0
    for i in 1:n
        min_i = minimum(cost_connection[i,:])
        if min_i > lb
            lb = min_i
        end
    end
    return lb
end

function one_centre_ub(n::Int, m::Int, cost_connection::Matrix{Int})
    """Compute an upper bound for 1 centre"""
    ub = maximum(cost_connection)
    for j in 1:m
        max_j = maximum(cost_connection[:,j])
        if max_j < ub
            ub = max_j
        end
    end
    return ub
end

function main_p_centres(n, m, cost_connection, nb_centres)
    # exemples d'utilisation de PLS
    # formulation forte, pb relache, pas d'affichage
    println("\n")
    println("Valeur de la borne inférieure = ", lower_bound_0(n, m, cost_connection), "\n")
    println("Valeur de la borne supérieure obtenue avec 1 centre = ", one_centre_ub(n, m, cost_connection), "\n")
    obj, temps, noeuds, x, y = solve_pc1(n, m, cost_connection, nb_centres, relaxation =  true, verbose = 0)
    println("Valeur relaxation formulation PC1 = ", obj, ", obtenue en ", temps, "\n")
    obj, temps, noeuds, x, y = solve_pc1(n, m, cost_connection, nb_centres, relaxation =  false, verbose = 0)
    println("Valeur formulation PC1 = ", obj,  ", obtenue en ", temps,"\n")
    obj, temps, noeuds, x, y = solve_pc(n, m, cost_connection, nb_centres, relaxation =  true, verbose = 0)
    println("Valeur relaxation formulation PC = ", obj, ", obtenue en ", temps, "\n")
    obj, temps, noeuds, x, y = solve_pc(n, m, cost_connection, nb_centres, relaxation =  false, verbose = 0)
    println("Valeur formulation PC = ", obj, ", obtenue en ", temps, "\n")
end