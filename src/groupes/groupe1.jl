# Valentin et Agathe

# Valentin et Agathe

using JuMP
using CPLEX
using DataFrames
using CSV
using DataFrames
using CSV
include("../readData.jl")


function PL(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int}; formulation::String = "S", pb_relache::Bool = false, silence::Bool = false)
    # n nombre de sites
    # m le nombre de clients
    #n, m = m, n
    # m nombre de sites
    # n le nombre de clients

    model = Model(CPLEX.Optimizer)
    if silence
        set_silent(model)
    end
    set_time_limit_sec(model, 1800)
    MOI.set(model, MOI.NumberOfThreads(), 1)

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

    sites = [i for i in 1:m if value(y[i]) > 1e-5]

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        if !pb_relache
            gap = MOI.get(model, MOI.RelativeGap())
        else 
            gap = 0
        end
        return objective_value(model), solve_time(model), node_count(model), value.(x), sites, gap
    end
    sites = [i for i in 1:m if value(y[i]) > 1e-5]

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        if !pb_relache
            gap = MOI.get(model, MOI.RelativeGap())
        else 
            gap = 0
        end
        return objective_value(model), solve_time(model), node_count(model), value.(x), sites, gap
    end
end


function main_PL(n, m, opening_cost, cost_connection)
    # exemples d'utilisation de PL
    # formulation forte, pb relache, pas d'affichage
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  true, silence = true)
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  true, silence = true)
    println("Valeur relaxation formulation forte = ", obj, "\n \n")

    # formulation faible, pb relache, pas d'affichage
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "W", pb_relache =  true, silence = true)
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "W", pb_relache =  true, silence = true)
    println("Valeur relaxation formulation faible = ", obj, "\n \n")

    # formulation forte, PLNE, pas d'affichage
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  false, silence = true)
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  false, silence = true)
    println("Valeur formulation forte = ", obj, "\n \n")
    
    # formulation faible, PLNE, pas d'affichage
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "W", pb_relache =  false, silence = true)
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "W", pb_relache =  false, silence = true)
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

    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  false, silence = true)
    println("cout PL = ", obj)
    println("capteurs PL = ", [i for (i, value) in enumerate(y) if value >= 1e-5])
    obj, temps, noeuds, x, y = PL(n, m, opening_cost, cost_connection, formulation =  "S", pb_relache =  false, silence = true)
    println("cout PL = ", obj)
    println("capteurs PL = ", [i for (i, value) in enumerate(y) if value >= 1e-5])

    #verification heurGlou
    #res = sum(opening_cost[j] for j in S)
    #for i in 1:n
    #    res += minimum([cost_connection[i,j] for j in S])
    #end
    #print("verification res = ", res)

end

function pipeline_heurGlou()
    folder_path = "./data/"
    output_file = "./src/resultats/groupe1/heurGlou.csv"

    if isfile(output_file) 
        println("fichier deja existant")
        #result_df = CSV.read(output_file, DataFrame, types=String)
        result_df = DataFrame(CSV.File(output_file))
        files = setdiff(readdir(folder_path), result_df[:, "instance"])
    else
        result_df = DataFrame(instance = [], UB=[], LB=[], temps=[], sites=[])
        files = readdir(folder_path)
    end

    for file in files 
        n, m, opening_cost, cost_connection = read_data(folder_path * file)

        start_time = time()
        sites_heurGlou, v_obj_heurGlou, nu = heurGlou(n, m, opening_cost, cost_connection)
        temps_heurGlou = time() - start_time

        push!(result_df, (instance = file, UB = v_obj_heurGlou, LB= sum(nu[i] for i in 1:n), temps=temps_heurGlou,sites=sites_heurGlou))       
        CSV.write(output_file, DataFrame(result_df))    
    end
end

function vec_to_string(vector::Vector{Int})
    str = join(string.(vector), ",")  # Convert each element of the vector to a string and join them with a comma and space separator
    return "[$str]"  # Surround the resulting string with square brackets
end

function pipeline_PL(formulation::String = "W")
    folder_path = "./data/" 
    output_file = "./results/PL" * formulation * ".csv"

    if isfile(output_file) 
        println("fichier deja existant")
        #result_df = CSV.read(output_file, DataFrame, types=String)
        result_df = DataFrame(CSV.File(output_file))
        result_df[!,:sites] = convert.( String255,result_df[!,:sites])
        files = setdiff(readdir(folder_path), result_df[:, "instance"])
        println("a faire tourner : ", files)
    else
        println("Pas de fichier existant")
        result_df = DataFrame(instance = [], obj=[], v_rel=[], gap = [], noeuds = [], temps=[], sites=[])
        files = readdir(folder_path)
    end
    #files = ["instRand_100_100_1.txt", "instTest.txt"]
    #files = []

    for file in files 
        n, m, opening_cost, cost_connection = read_data(folder_path * file)

        v_rel, _, _, _, _, _ = PL(n, m, opening_cost, cost_connection, formulation =  formulation, pb_relache =  true, silence = false)   
        obj, temps, noeuds, _, sites, gap = PL(n, m, opening_cost, cost_connection, formulation =  formulation, pb_relache =  false, silence = false)

        push!(result_df, (instance = file, obj = obj, v_rel = v_rel, gap = gap, noeuds = noeuds, temps = temps, sites =vec_to_string(sites)))       
        CSV.write(output_file, DataFrame(result_df))    
    end
end

function main_pipeline()
    pipeline_PL("S")
    pipeline_PL("W")
end