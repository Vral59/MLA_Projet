# Justin et Benoit
# 2e formulation MILP sur le problème de Localisation simple

using JuMP
using CPLEX
using CSV
using DataFrames
include("../readData.jl")


function benchmark_grp2()
    repo_path = "data"
    time_limit = 60
    silence = true

    columns = ["instance", "method", "obj", "bound", "time", "n_nodes",  "root_obj", "root_bound", "n_variables", "n_constraints"]
    data = []
    for entry in readdir(repo_path)
        fullpath = joinpath(repo_path, entry)
        if isfile(fullpath)
            println("Résolution pour le fichier : $fullpath")
        end
        n, m, opening_cost, cost_connection = read_data(fullpath)

        # Formulation alternative
        println("\nFormulation alternative relachée")
        root_obj, root_bound, _, _, _, _ = PLS_bis(n, m, opening_cost,
        cost_connection, pb_relache = true, silence = silence, time_limit = time_limit)

        println("\nFormulation alternative exacte")
        obj, bound, time, n_nodes, n_variables, n_constraints = PLS_bis(n, m, opening_cost,
        cost_connection, pb_relache = false, silence = silence, time_limit = time_limit)
        push!(data, (entry, "formul_alt", obj, bound, time, n_nodes, root_obj, root_bound, n_variables, n_constraints))

        # Formulation alternative variante
        println("\nFormulation alternative variante relachée")
        root_obj, root_bound, _, _, _, _ = PLS_bis(n, m, opening_cost,
        cost_connection, pb_relache = true, variante = true, silence = silence, time_limit = time_limit)

        println("\nFormulation alternative variante exacte")
        obj, bound, time, n_nodes, n_variables, n_constraints = PLS_bis(n, m, opening_cost,
            cost_connection, pb_relache = false, variante = true, silence = silence, time_limit = time_limit)
        push!(data, (entry, "formul_alt_variante", obj, bound, time, n_nodes, root_obj, root_bound, n_variables, n_constraints))

        df = DataFrame(data, columns)
        CSV.write("results/benchmark_grp2_.csv", df)
    end

end


function main_pls_bis(n, m, opening_cost, cost_connection; time_limit = 30, silence = true)
    # Formulation alternative, relaché
    obj, bound, time, n_nodes, n_variables, n_constraints = PLS_bis(n, m, opening_cost,
        cost_connection, pb_relache = true, silence = silence, time_limit = time_limit)
    println("\nValeur relaxation formulation alternative = ", round(obj, digits=1), " en " , round(time, digits=3), "s. ",
        "Noeuds explorés: ", n_nodes, ". Meilleure borne: ", round(bound, digits=1))

    # Formulation alternative, exacte
    obj, bound, time, n_nodes, n_variables, n_constraints = PLS_bis(n, m, opening_cost,
        cost_connection, pb_relache = false, silence = silence, time_limit = time_limit)
    println("Valeur formulation alternative = ", round(obj, digits=1), " en " , round(time, digits=3), "s. ",
        "Noeuds explorés: ", n_nodes, ". Meilleure borne: ", round(bound, digits=1))

    # Formulation alternative, relaché, variante
    obj, bound, time, n_nodes, n_variables, n_constraints = PLS_bis(n, m, opening_cost,
        cost_connection, pb_relache = true, variante = true, silence = silence, time_limit = time_limit)
    println("Valeur relaxation formulation alternative (variante) = ", round(obj, digits=1), " en " , round(time, digits=3), "s. ",
        "Noeuds explorés: ", n_nodes, ". Meilleure borne: ", round(bound, digits=1))

    # Formulation alternative, exacte, variante
    obj, bound, time, n_nodes, n_variables, n_constraints = PLS_bis(n, m, opening_cost,
        cost_connection, pb_relache = false, variante = true, silence = silence, time_limit = time_limit)
    println("Valeur formulation alternative (variante) = ", round(obj, digits=1), " en " , round(time, digits=3), "s. ",
        "Noeuds explorés: ", n_nodes, ". Meilleure borne: ", round(bound, digits=1))
end


function PLS_bis(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int};
        pb_relache::Bool = false, silence::Bool = true, variante::Bool = false, time_limit::Int = 0)
    """2e formulation MILP sur le problème de Localisation simple

    - m nombre de sites
    - n le nombre de clients

    Pour chaque site j:
    y[j] = 1 si le site j est ouvert, 0 sinon

    Pour chaque client i:
    - on note k_i le nombre de distances distinctes entre le client i et les sites
    - on crée également des C_i tels que C_i^1 < C_i^2 < ... < C_i^{k_i}

    A la place des variables x[i,j] du problème initial, on introduit une variable z_i
    vecteur de booléens pour chaque client i, de taille k_i.
    z[i,k] = 1 si aucun site n'est ouvert à une distance <= C_i^k du client i

    z[i,k_i] = 0. Tous les clients sont servis par au moins un site
    """
    if !silence
        println("Building problem...")
    end

    model = Model(CPLEX.Optimizer)
    MOI.set(model, MOI.RelativeGapTolerance(), 1e-6)

    if silence
        set_silent(model)
    end
    if time_limit > 0
        set_time_limit_sec(model, time_limit)
    end

    sorted_distances = map(i -> sort(unique(cost_connection[i, :])), 1:n)
    max_distinct_distances = maximum(length.(sorted_distances))
    # map of vector. If the problem is large, would be better to be an array and padd with anything at the end

    @variable(model, y[1:m], Bin)
    @variable(model, z[1:n, 1:max_distinct_distances], Bin)

    @constraint(model, [i in 1:n], z[i,1] + sum(y.*(cost_connection[i,:] .== sorted_distances[i][1])) >= 1)
    if variante
        @constraint(model, [i in 1:n, k in 2:length(sorted_distances[i])],
            z[i,k] - z[i,k-1] + sum(y.*(cost_connection[i,:] .== sorted_distances[i][k])) >= 0)
    else
        @constraint(model, [i in 1:n, k in 2:length(sorted_distances[i])],
            z[i,k] + sum(y.*(cost_connection[i,:] .<= sorted_distances[i][k])) >= 1)
    end
    @constraint(model, [i in 1:n], z[i, length(sorted_distances[i]):end] .== 0)


    opening_cost_value = sum(opening_cost.*y)
    # Large problems are long to write, because there are too many "+" operations in the following lines (for)
    # it would be better if sorted_distances was an array so that operations between arrays can be used
    assignation_cost = begin
        sum(sorted_distances[i][1] for i in 1:n) +
        sum(z[i,k]*(sorted_distances[i][k+1] - sorted_distances[i][k]) for i in 1:n for k in 1:length(sorted_distances[i])-1)
    end
    total_cost = opening_cost_value + assignation_cost
    @objective(model, Min, total_cost)

    if pb_relache
        relax_integrality(model)
    end

    if !silence
        println("Problem built! Starting optimization...")
    end

    optimize!(model)

    n_variables = num_variables(model)
    n_constraints = sum(num_constraints(model, F, S) for (F, S) in list_of_constraint_types(model))
    n_nodes = node_count(model)
    if has_values(model)
        obj_value = objective_value(model)
    else
        obj_value = 1e8
    end
    lower_bound = objective_bound(model)
    time = solve_time(model)

    if !silence
        println("Nombre de variables: ", n_variables)
        println("Nombre de contraintes: ", n_constraints)

        println("Objective value: ", obj_value)
        if has_values(model)
            println("Opening cost: ", value(opening_cost_value))
            println("Assignation cost: ", value(assignation_cost))
        end

        println("Lower bound: ", lower_bound)
        println("Solve time: ", time)
        println("Number of nodes: ", n_nodes)
    end

    if pb_relache && n_nodes != 0
        println("Erreur: La relaxation doit être résolue en un noeud")
    end

    return obj_value, lower_bound, time, n_nodes, n_variables, n_constraints
end
