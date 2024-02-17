using JuMP
using CPLEX
include("readData.jl")


function PLS(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int}; formulation::String = "S", pb_relache::Bool = false, silence::Bool = false)
    """formulation forte"""
    # n nombre de sites
    # m le nombre de clients
    #n, m = m, n
    # m nombre de sites
    # n le nombre de clients

    model = Model(CPLEX.Optimizer)
    if silence
        set_silent(model)
    end

    @variable(model, x[1:n, 1:m] >= 0)
    @variable(model, y[1:m], Bin)

    # choix de la formulation (S pour fort, W pour faible)
    if formulation == "S"
        println("formulation forte")
        @constraint(model, [i in 1:n, j in 1:m], x[i,j] <= y[j])
    else 
        if formulation == "W"
            println("formulation faible")
            @constraint(model, [j in 1:m], sum(x[i,j] for i in 1:n) <= y[j] * n )
        else
            @warn("Mauvais appel de la formulation : mettre 'S' pour la formulation forte, 'W' pour la formulation faible")
            return(-1)
        end
    end

    @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:m)  == 1)

    # cost connection s'appelle de manière différente de l'enonce
    @objective(model, Min, sum(cost_connection[i,j] * x[i,j] for i in 1:n for j in 1:m) + sum(opening_cost[j] * y[j] for j in 1:m))

    if pb_relache
        relax_integrality(model)
    end

    optimize!(model)

    return objective_value(model), solve_time(model), node_count(model), value.(x), value.(y)
end


function main()
    n, m, opening_cost, cost_connection = read_data("data/instTest.txt")

    # exemples d'utilisation de PLS
    # formulation forte, pb relache, pas d'affichage
    obj, temps, noeuds, x, y = PLS(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  true, silence = true)
    println("Valeur relaxation formulation forte = ", obj, "\n \n")

    # formulation faible, pb relache, pas d'affichage
    obj, temps, noeuds, x, y = PLS(n, m, opening_cost, cost_connection, formulation =  "W", pb_relache =  true, silence = true)
    println("Valeur relaxation formulation faible = ", obj, "\n \n")

    # formulation forte, PLNE, pas d'affichage
    obj, temps, noeuds, x, y = PLS(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  false, silence = true)
    println("Valeur formulation forte = ", obj, "\n \n")
    
    # formulation faible, PLNE, pas d'affichage
    obj, temps, noeuds, x, y = PLS(n, m, opening_cost, cost_connection, formulation =  "W", pb_relache =  false, silence = true)
    println("Valeur formulation faible = ", obj)

    # les trois derniers arguments peuvent ne pas être specifiés. valeurs par défaut: S, false, false (formulation forte, PLNE, modele bavard)
    
   
end