# Groupe Francis et Jeanne

using JuMP
using CPLEX
# using CSV
using DataFrames

include("groupe4.jl")


function NPCiR(n::Int, m::Int, p::Int, distances::Matrix{Int}, nb_z_fixed = Int)
    """formulation NPCi relâché"""
    # y relâché

    # n nombre de clients
    # m nombre de sites
    # p nombre de sites à ouvrir
    # distances matrice des distances/coûts

    distances_K_D = distances_triées(n::Int, m::Int, distances)
    K, D = distances_K_D[1], distances_K_D[2]


    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model,"CPX_PARAM_SCRIND",0)

    @variable(model, 0<=z[1:K]<=1)
    @variable(model, 0<=y[1:m]<=1)

    # Les nb_z_fixed premiers z sont fixés à 1
    @constraint(model,[i in 1:nb_z_fixed],z[i]==1)

    
    @constraint(model, [i in 1:n], sum(y[j] for j in 1:m)  == p)
    
    @constraint(model, [i in 1:n, k in 1:K], z[k]+sum(y[j] for j in 1:m if distances[i,j]<=D[k])  >= 1)

    @constraint(model, [k in 1:(K-1)], z[k] >= z[k+1])



    @objective(model, Min, D[1]+sum((D[k]-D[k-1]) * z[k-1] for k in 2:K))

    optimize!(model)


    return value.(z)
end

function algo_fixation_z(n::Int, m::Int, p::Int, distances::Matrix{Int})
    """Algorithme de fixation des z"""
    stop = false
    nb_z_fixed = 0 #init

    while stop == false
        z = NPCiR(n, m, p, distances, nb_z_fixed)

        add_fixed_z = 0 # nombre de z à fixer supplémentaires

        for z_i in value.(z)[nb_z_fixed+1:end]
            if z_i==1
                add_fixed_z += 1
            end               
        end

        if value.(z)[nb_z_fixed + add_fixed_z + 1] == 0
            stop = true
        else
            add_fixed_z += 1
        end

        nb_z_fixed += add_fixed_z
        println("nb_z_fixed = ", nb_z_fixed)
    end

    return distances_triées(n, m, distances)[2][nb_z_fixed]
end

