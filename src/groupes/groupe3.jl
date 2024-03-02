# Guillaume Casanova et Raphael Taisant
# P-centres problems


function solve_pc1(n::Int, m::Int, cost_connection::Matrix{Int}, nb_centres::Int; relaxation::Bool=false, time_limit::Float64=0, verbose::Int=0)
    """Solve the p-centers problem with formulation PC1"""
    # n number of clients
    # m number of sites

    model = Model(CPLEX.Optimizer)
    if verbose == 0
        set_silent(model)
    end
    if time_limit > 0
        set_time_limit_sec(model, time_limit)
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

    n_variables = num_variables(model)
    n_constraints = sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model))

    return objective_value(model), solve_time(model), node_count(model), value.(x), value.(y), n_variables, n_constraints, objective_bound(model)
end

function solve_pc(n::Int, m::Int, cost_connection::Matrix{Int}, nb_centres::Int; relaxation::Bool=false,  time_limit::Float64=0, verbose::Int=0)
    """Solve the p-centers problem with formulation PC1"""
    # n number of clients
    # m number of sites

    model = Model(CPLEX.Optimizer)
    if verbose == 0
        set_silent(model)
    end
    if time_limit > 0
        set_time_limit_sec(model, time_limit)
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

    n_variables = num_variables(model)
    n_constraints = sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model))

    return objective_value(model), solve_time(model), node_count(model), value.(x), value.(y), n_variables, n_constraints, objective_bound(model)
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


function open_best_center!(x::Matrix{Int64}, y::Vector{Int64}, d_min::Union{Vector{Int64}, Vector{Float64}}, n::Int, m::Int, cost_connection::Matrix{Int})
    """Open the best site, modify the current solution x, y, d_min inplace.
    d_min is the distance to the closest center for each client"""
    best_obj, best_j = Inf, nothing
    clients_served = [[] for j in 1:m] # for a given facility j, used to update d_min without copying the whole array
    for j in 1:m
        if y[j] == 1
            continue
        end
        obj = 0
        for i in 1:n
            d = min(cost_connection[i,j], d_min[i])
            obj = max(obj, d)
            if cost_connection[i,j] < d_min[i]
                push!(clients_served[j], i)
            end
        end
        if obj < best_obj
            best_obj = obj
            best_j = j
        end
    end
    
    # Modify the solution
    y[best_j] = 1
    for i in clients_served[best_j]
        d_min[i] = cost_connection[i,best_j]
        x[i,best_j] = 1
    end

    @assert best_obj == maximum(d_min)
    return best_obj, best_j
end

function init_empty_solution(n::Int, m::Int)
    """Return an empty solution"""
    x = zeros(Int, n, m)
    y = zeros(Int, m)
    d_min = [Inf for _ in 1:n]
    return x, y, d_min
end

function  solve_greedy(n::Int, m::Int, cost_connection::Matrix{Int}, nb_centres::Int; verbose::Int = 0)
    """Compute a solution with the greedy heuristic adding centers sequentially"""
    obj = Inf
    x, y, d_min = init_empty_solution(n, m)
    for _ in 1:nb_centres
        obj, j = open_best_center!(x, y, d_min, n, m, cost_connection)
        if verbose == 1 println("Ouverture du site $j, nouvelle valeur = $obj") end
    end
    return obj, x, y
end

function main_p_centres(n::Int, m::Int, cost_connection::Matrix{Int64}, nb_centres::Int)
    # exemples d'utilisation de PLS
    # formulation forte, pb relache, pas d'affichage
    println("\n")
    println("Valeur de la borne inférieure = ", lower_bound_0(n, m, cost_connection), "\n")
    println("Valeur de la borne supérieure obtenue avec 1 centre = ", one_centre_ub(n, m, cost_connection), "\n")
    obj, temps, noeuds, x, y, _, _, _ = solve_pc1(n, m, cost_connection, nb_centres, relaxation =  true, time_limit=60, verbose = 0)
    println("Valeur relaxation formulation PC1 = ", obj, ", obtenue en ", temps, "s", "\n")
    obj, temps, noeuds, x, y, _,_,_ = solve_pc1(n, m, cost_connection, nb_centres, relaxation =  false, time_limit=60, verbose = 0)
    println("Valeur formulation PC1 = ", obj,  ", obtenue en ", temps, "s","\n")
    obj, temps, noeuds, x, y,_,_,_ = solve_pc(n, m, cost_connection, nb_centres, relaxation =  true, time_limit=60, verbose = 0)
    println("Valeur relaxation formulation PC = ", obj, ", obtenue en ", temps, "s", "\n")
    obj, temps, noeuds, x, y,_,_,_ = solve_pc(n, m, cost_connection, nb_centres, relaxation =  false, time_limit=60, verbose = 0)
    println("Valeur formulation PC = ", obj, ", obtenue en ", temps, "s", "\n")
    temps = @elapsed obj, x, y = solve_greedy(n, m, cost_connection, nb_centres, verbose = 1)
    println("Valeur heuristique greedy = ", obj, ", obtenue en ", temps, "s", "\n")
end

function benchmark_grp3()
    repo_path = "tsp_data"
    time_lim = 60.0

    ps = [5, 10, 20, 50]

    columns = ["instance", "p", "method", "obj", "bound", "time", "n_nodes", "root_obj", "n_variables", "n_constraints"]
    data = []
    for entry in readdir(repo_path)
        # entry = "dj38.tsp"
        fullpath = joinpath(repo_path, entry)
        if isfile(fullpath)
            println("Résolution pour le fichier : $fullpath")
        end
        n, m, distances = readInstance_tsp(fullpath)

        for p in ps
            if p > n
                break
            end

            # PC1
            println("PC1 relax")
            root_obj, _, _, _, _, _, _,_ = solve_pc1(n, m, distances, p, relaxation = true, time_limit=time_lim, verbose = 0)
            println("PC1")
            obj, time, n_nodes, _, _, n_variables, n_constraints, lb = solve_pc1(n, m, distances, p, relaxation = false, time_limit=time_lim, verbose = 0)
            push!(data, (entry, p, "PC1", round(obj, digits=1), round(lb, digits=1), round(time, digits=3),
                n_nodes, round(root_obj, digits=1), n_variables, n_constraints))
            
            # PC
            println("PC relax")
            root_obj, _, _, _, _, _, _,_ = solve_pc(n, m, distances, p, relaxation = true, time_limit=time_lim, verbose = 0)
            println("PC")
            obj, time, n_nodes, _, _, n_variables, n_constraints, lb = solve_pc(n, m, distances, p, relaxation = false, time_limit=time_lim, verbose = 0)
            push!(data, (entry, p, "PC", round(obj, digits=1), round(lb, digits=1), round(time, digits=3),
                n_nodes, round(root_obj, digits=1), n_variables, n_constraints))
            
            # Greedy
            println("Greedy")
            time = @elapsed obj, _, _,= solve_greedy(n, m, distances, p, verbose = 0)
            push!(data, (entry, p, "Greedy", round(obj, digits=1), round(obj, digits=1), round(time, digits=3),
                0, 0, 0, 0))

            df = DataFrame(data, columns)
            CSV.write("results/benchmark_grp3_.csv", df)
        end
        break
    end
end