# Justin et Benoit
# 2e formulation MILP sur le problème de Localisation simple

using JuMP
using CPLEX
using LinearAlgebra
include("../readData.jl")

function PLS_bis(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int};
        formulation::String = "S", pb_relache::Bool = false, silence::Bool = false)
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

    @variable(model, y[1:n], Bin)

    sorted_distances = []
    max_distinct_distances = 0
    for i in 1:n
        dist_i = sort(unique(cost_connection[i,:]))
        sorted_dist_i = sort(dist_i)
        k_i = length(dist_i)

        if k_i > max_distinct_distances
            max_distinct_distances = k_i
        end
        push!(sorted_distances, sorted_dist_i)
    end

    @variable(model, z[1:n, 1:max_distinct_distances], Bin)

    ## To complete
end