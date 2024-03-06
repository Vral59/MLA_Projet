# Tiphaine George et Louise Lallemand

using JuMP
using CPLEX
using DataFrames
using CSV
include("../readData.jl")

function PLR(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int}; silence::Bool = false)
    # n nombre de sites
    # m le nombre de clients
    #n, m = m, n
    # m nombre de sites
    # n le nombre de clients

    model = Model(CPLEX.Optimizer)
    if silence
        set_silent(model)
    end
    set_time_limit_sec(model, 600)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    @variable(model, x[1:n, 1:m] >= 0)
    @variable(model, 0<= y[1:m] <=1)

    # contrainte forte
    cons_w = @constraint(model, [i in 1:n, j in 1:m], x[i,j] <= y[j])

    cons_v = @constraint(model, [i in 1:n], sum(x[i,j] for j in 1:m)  == 1)

    # cost connection s'appelle de manière différente de l'enonce
    @objective(model, Min, sum(cost_connection[i,j] * x[i,j] for i in 1:n for j in 1:m) + sum(opening_cost[j] * y[j] for j in 1:m))
    optimize!(model)

    sites = [i for i in 1:m if value(y[i]) > 1e-5]
    constraints = [constr for constr in all_constraints(model; include_variable_in_set_constraints=false)]
    dual_values = [dual(constr) for constr in constraints]
    dual_v = [dual(constr) for constr in cons_v]
    dual_w = [dual(constr) for constr in cons_w]

    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        return objective_value(model), solve_time(model), value.(x), sites, dual_v, dual_w
    end
end


function methode_arrondi(n, m, opening_cost, cost_connection, x, y, v, w)
    # à partir des résultats du PLR, construit les clusters puis leur affecte un site 
    t0 = time()
    G = build_G(x)
    Clusters = partition(G,v)
    xh,yh = affectation(Clusters, G, opening_cost)
    valh = objectif_value(xh,yh,opening_cost,cost_connection)
    t1 = time()
    return valh, sites(yh), t1-t0
end

function sites(y)
    m = size(y,1)
    sites = Vector{Int}()
    for s in 1:m
        if y[s]==1
            push!(sites, s)
        end
    end
    return sites
end

function build_G(x)
    # G = matrice d'adjacence, Gij=1 ssi xij>0
    n = size(x,1)
    m = size(x,2)
    G = zeros(n,m)
    for i in 1:n
        for j in 1:m
            if x[i,j] > 0
                G[i,j] = 1
            end
        end
    end
    return G
end

function partition(G, v)
    vmax = maximum(v)
    n = size(G,1)
    Clusters = []
    affecte = [false for _ in 1:n]
    while (false in affecte)
        cand = vmax+1
        i_cand = 0
        for i in 1:n
            if !affecte[i] && v[i]<cand
                cand = v[i]
                i_cand = i
            end
        end
        voisins = find_voisins(G,i_cand,affecte)
        for vois in voisins
            affecte[vois] = true
        end
        push!(Clusters, voisins)
    end
    return Clusters
end

function find_voisins(G, i, affecte)
    n = size(G,1)
    m = size(G,2)
    voisins = [i]
    for j in 1:m
        if G[i,j]==1
            for k in 1:n
                if !affecte[k] && G[k,j]==1
                    push!(voisins, k) 
                end
            end
        end
    end
    return unique(voisins)
end

function affectation(Clusters, G, opening_cost)
    n = size(G,1)
    m = size(G,2)
    fmax = maximum(opening_cost)
    x = zeros(n,m)
    y = zeros(m)
    for cl in Clusters
        ik = cl[1]
        j_cand = 0
        f_cand = fmax+1 
        for j in 1:m
            if G[ik,j]==1 && opening_cost[j]<f_cand
                f_cand = opening_cost[j]
                j_cand = j
            end
        end
        y[j_cand] = 1
        for i in cl
            x[i,j_cand] = 1
        end
    end
    return x,y
end

function objectif_value(x,y,opening_cost,cost_connection)
    n = size(x,1)
    m = size(x,2)
    v=sum(cost_connection[i,j]*x[i,j] for i in 1:n for j in 1:m) + sum(opening_cost[j]*y[j] for j in 1:m)
    return v
end

# adapté de Agathe et Valentin
function pipeline_MA()
    folder_path = "./data/"
    output_file = "./src/resultats/groupe7/methode_arrondi.csv"

    if isfile(output_file) 
        println("fichier deja existant")
        #result_df = CSV.read(output_file, DataFrame, types=String)
        result_df = DataFrame(CSV.File(output_file))
        result_df[!,:sites] = convert.( String255,result_df[!,:sites])
        files = setdiff(readdir(folder_path), result_df[:, "instance"])
    else
        result_df = DataFrame(instance = [], obj=[], v_rel=[], temps=[], sites=[])
        files = readdir(folder_path)
    end
    #files = ["instRand_100_100_1.txt", "instTest.txt"]
    #files = []

    for file in files 
        n, m, opening_cost, cost_connection = read_data(folder_path * file)

        obj, temps_pl, x, y, v, w = PLR(n, m, opening_cost, cost_connection, silence = true)
        valh, sites, temps_ma = methode_arrondi(n, m, opening_cost, cost_connection, x, y, v, w)
        push!(result_df, (instance = file, obj = valh, v_rel = obj, temps = temps_ma+temps_pl, sites =vec_to_string(sites)))       
        CSV.write(output_file, DataFrame(result_df))    
    end
end

# copié de Agathe et Valentin
function vec_to_string(vector::Vector{Int})
    str = join(string.(vector), ",")  # Convert each element of the vector to a string and join them with a comma and space separator
    return "[$str]"  # Surround the resulting string with square brackets
end