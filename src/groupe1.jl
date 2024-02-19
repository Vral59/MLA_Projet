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


function main_pls(n, m, opening_cost, cost_connection)
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

function partiePos(a)
    if a >= 0
        return a
    end
    return 0
end

function heurGlou(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int})
    """heurGlou retourne S, z, nu comme définis dans l'algorithme Chapitre 1"""
    j_min = argmin([opening_cost[j] + sum(cost_connection[i,j] for i in 1:n) for j in 1:m])
    S = [j_min]
    nu = [cost_connection[i,j_min] for i in 1:n]
    z = opening_cost[j_min] + sum(nu[i] for i in 1:n)
    stop = false
    while !stop
        println("S = ", S)
        sites = setdiff(1:m, S)
        Deltas = [opening_cost[j] - sum(partiePos(nu[i] - cost_connection[i,j]) for i in 1:n) for j in sites]
        j_min = argmin(Deltas) 
        site_j_min = sites[argmin(Deltas)] # sinon probleme d'indice comme Delta n'a pas tous les indices
        if Deltas[j_min] >= 0
            return(S, z, nu)
        else
            z += Deltas[j_min]
            push!(S, site_j_min)
            nu = [nu[i] - partiePos(nu[i] - cost_connection[i,site_j_min]) for i in 1:n]
        end   
    end
end


function main_heurGlou(n, m, opening_cost, cost_connection)
    """heurGlou retourne S, z, nu comme définis dans l'algorithme Chapitre 1"""
    S, z, nu = heurGlou(n, m, opening_cost, cost_connection)
    println("cout heurGlou = ", z)
    println("capteurs heurGlou = ", S, "\n")

    obj, temps, noeuds, x, y = PLS(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  false, silence = true)
    println("cout PLS = ", obj)
    println("capteurs PLS = ", [i for (i, value) in enumerate(y) if value >= 1e-5])

    #verification heurGlou
    #res = sum(opening_cost[j] for j in S)
    #for i in 1:n
    #    res += minimum([cost_connection[i,j] for j in S])
    #end
    #print("verification res = ", res)

end
