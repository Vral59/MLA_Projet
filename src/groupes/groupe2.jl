# Justin et Benoit
# 2e formulation MILP sur le problème de Localisation simple

using JuMP
using CPLEX
using LinearAlgebra
include("../readData.jl")

function main_pls_bis(n, m, opening_cost, cost_connection)
    # Formulation alternative, relaché
    obj, time, node, z, y = PLS_bis(n, m, opening_cost, cost_connection, pb_relache = true)
    println("\nValeur relaxation formulation alternative = ", obj, "\n \n")
    # Formulation alternative, exacte
    obj, time, node, z, y = PLS_bis(n, m, opening_cost, cost_connection, pb_relache = false)
    println("Valeur formulation alternative = ", obj, "\n \n")
    # Formulation alternative, relaché, variante
    obj, time, node, z, y = PLS_bis(n, m, opening_cost, cost_connection, pb_relache = true, variante = true)
    println("Valeur relaxation formulation alternative (variante) = ", obj, "\n \n")
    # Formulation alternative, exacte, variante
    obj, time, node, z, y = PLS_bis(n, m, opening_cost, cost_connection, pb_relache = false, variante = true)
    println("Valeur formulation alternative (variante) = ", obj, "\n \n")
end

function PLS_bis(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int}; pb_relache::Bool = false, silence::Bool = true, variante::Bool = false)
    """2e formulation MILP sur le problème de Localisation simple

    - m nombre de sites
    - n le nombre de clients

    Pour chaque client i:
    - on note k_i le nombre de distances distinctes entre le client i et les sites
    - on crée également des C_i tels que C_i^1 < C_i^2 < ... < C_i^{k_i}

    A la place des variables x[i,j] du problème initial, on introduit une variable z_i
    vecteur de booléens pour chaque client i, de taille k_i.
    z[i,k] = 1 si aucun site n'est ouvert à une distance <= C_i^k du client i
    """

    model = Model(CPLEX.Optimizer)
    if silence
        set_silent(model)
    end

    @variable(model, y[1:m], Bin)

    sorted_distances = []
    max_distinct_distances = 0
    for i in 1:n
        sorted_dist_i = sort(unique(cost_connection[i,:]))
        k_i = length(sorted_dist_i)

        if k_i > max_distinct_distances
            max_distinct_distances = k_i
        end
        push!(sorted_distances, sorted_dist_i)
    end

    @variable(model, z[1:n, 1:max_distinct_distances], Bin)

    ## To complete

    for i in 1:n
        @constraint(model, z[i,1] + sum(y[j] for j in 1:m if cost_connection[i,j] == sorted_distances[i][1]) >= 1)
    end

    if variante
        for i in 1:n
            @constraint(model, [k in 2:length(sorted_distances[i])], z[i,k] + z[i,k-1] + sum(y[j] for j in 1:m if cost_connection[i,j] == sorted_distances[i][k]) >= 1)
        end
    else
        for i in 1:n
            @constraint(model, [k in 2:length(sorted_distances[i])], z[i,k] + sum(y[j] for j in 1:m if cost_connection[i,j] <= sorted_distances[i][k]) >= 1)
        end
    end

    @constraint(model, [i in 1:n, k in 1:max_distinct_distances; k >= length(sorted_distances[i])], z[i, k] == 0)

    @objective(model, Min, sum(opening_cost[j]*y[j] for j in 1:m) + sum(sorted_distances[i][1] for i in 1:n)
    +sum(z[i,k]*(sorted_distances[i][k+1]- sorted_distances[i][k]) for i in 1:n, k in 1:length(sorted_distances[i])-1))

    if pb_relache
        relax_integrality(model)
    end

    optimize!(model)

    return objective_value(model), solve_time(model), node_count(model), value.(z), value.(y)
end