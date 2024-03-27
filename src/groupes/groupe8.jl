# Ines Chergui et Margot Boyer
using JuMP
using CPLEX
using Dates
using Statistics
using StatsBase

include("../readData.jl")

current_dir = pwd()

path = joinpath(current_dir, "tsp_data\\uy734.tsp")
#path = "MLA_Projet\\data\\ga250a-4.txt"


# n,m,distances = readInstance_tsp(path)




function pc(p::Int,m::Int,I::Vector{Int},distances::Matrix{Int})
    """
    Calcul du PC avec le sous-ensemble de client I
    """
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    @variable(model, x[i in I, 1:m], Bin)
    @variable(model, y[1:m], Bin)
    @variable(model, r>=0)

    @constraint(model, [i in I, j in 1:m], x[i,j] <= y[j])
    @constraint(model, [i in I], sum(x[i,j] for j in 1:m)  == 1)
    @constraint(model, sum(y) == p)
    @constraint(model, [i in I, j in 1:m], distances[i,j] * x[i,j] <= r)
    @objective(model, Min, r)
    optimize!(model)

    n_variables = num_variables(model)
    n_constraints = num_constraints(model; count_variable_in_set_constraints=false) # Count only structural constraints

    status = termination_status(model)
    @assert status == MOI.OPTIMAL || status == MOI.TIME_LIMIT
    
    obj, x_val, y_val = missing, missing, missing
    if primal_status(model) == MOI.FEASIBLE_POINT
        obj = objective_value(model)
        x_val, y_val = value.(x), value.(y)
    end
    lb = missing
    try lb = objective_bound(model) catch end
    return obj, solve_time(model), node_count(model), x_val, y_val, n_variables, n_constraints, lb
end



function rayon(n::Int,m::Int,distances::Matrix{Int},y::Vector{Float64})
    """
    Calcul du rayon optimal lorsque les usines sont ouvertes selon y
    """
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    @variable(model, x[1:n, 1:m] >= 0, Bin)
    @variable(model, r>=0)

    @constraint(model, [i in 1:n, j in 1:m], x[i,j] <= y[j])
    @constraint(model, [i in 1:n], sum(x[i,:])  == 1)

    @constraint(model, [i in 1:n], sum(distances[i,:] .* x[i,:]) <= r)
    @objective(model, Min, r)

    optimize!(model)

    status = termination_status(model)
    @assert status == MOI.OPTIMAL || status == MOI.TIME_LIMIT

    obj, x_val = missing, missing
    if primal_status(model) == MOI.FEASIBLE_POINT
        obj = objective_value(model)
        x_val = value.(x)
    end
    return obj, x_val
end

function mise_a_jour(n::Int, m::Int, I::Vector{Int}, x_val::Matrix{Float64}, distances::Matrix{Int}, k::Int = 1, method::String = "Distance_max")
    """
    Met à jour le sous-ensemble de clients selon différentes méthodes
    Args : 
        - k : nombre de clients à ajouter au sous-ensemble I
    """
    Ibar = collect(setdiff(1:n, I))
    if size(Ibar)[1] < k
        clients = Ibar
    elseif method == "Distance_max"
        """
        Ajoute les k clients qui maximisent la distance aux usines ouvertes courantes
        """
        r = x_val .* distances
        r[I,:] .= -1
        distances_max = maximum(r, dims =  2)
        clients = sortperm(vec(distances_max), rev=true)[1:k]

    elseif method == "Aleatoire"
        """
        Ajoute k clients aléatoirement
        """
        clients = sample(Ibar,k,replace = false)
    end
    append!(I,clients)
    return I
end


function initialise_I(n::Int,m::Int,distances::Matrix{Int},k::Int=1,method="Aleatoire")
    """ Initialise le sous-ensemble de clients I
    Args : 
        - k : nombre de clients de l'ensemble d'inialisation I
    """
    if method == "Distance_moyenne"
        """
        Initialise I par les k premiers clients avec des distances aux usines maximales en moyenne 
        """
        distances_moyenne = mean(distances, dims =  2)
        clients_tries = sortperm(vec(distances_moyenne), rev=true)
        I = Vector{Int}(clients_tries[1:k])
    elseif method == "Distance_maximale"
        """
        Initialise I par les k premiers clients avec les plus grandes distances maximales aux usines
        """
        distances_max = maximum(distances, dims =  2)
        clients_tries = sortperm(vec(distances_max), rev=true)
        I = Vector{Int}(clients_tries[1:k])
    elseif method == "Aleatoire"
        """
        Initialise I avec k clients aléatoires
        """
        I = sample(1:n,k,replace = false)
    end
    return I
end



function algo(n::Int, m::Int, distances::Matrix{Int}, p::Int, k::Int=1, methods::Vector{String} = ["Distance_moyenne","Distance_max"])
    start_time = now()
    I = initialise_I(n,m,distances,k, methods[1])
    println("methods : ", methods)
    println("I : ", I)
    lb = 0
    ub = Inf
    while ub > lb
        vopt, time, node_count, x_val, y_val, n_variables, n_constraints, lb = pc(p,m,I,distances)
        lb = vopt
        println("Temps passe solve pc : ",time )
        ry, x_val = rayon(n,m,distances,y_val)
        ub = min(ub, ry)
        I = mise_a_jour(n, m, I, x_val, distances, k, methods[2])
        #println("Nouvel ensemble I : ",I)
        #println("Nouvelle borne inf : ", lb)
        #println("Nouvelle borne sup : ", ub)
    end
    execution_time = Dates.value(now() - start_time)/1000
    println("Temps d'exécution(s) : ", execution_time)
    return lb, execution_time
end



function benchmark_grp8()

    current_dir = pwd()
    repo_path = "tsp_data\\"
    
    test_files = Dict("dj38.tsp" => [5,10,20], "wi29.tsp" => [5,10,20], "qa194.tsp" => [5,10])
    k_list = [1,3,5]

    methods_init = ["Aleatoire","Distance_maximale","Distance_moyenne"]
    methods_mise_a_jour = ["Aleatoire","Distance_max"]

    out = open("results/benchmark_grp8.csv","w")
    write(out,"instance;n;m;p;k;initialisation;mise_a_jour;obj;time\n")

    for instance in ["dj38.tsp","wi29.tsp","qa194.tsp"]
        println("Instance : ", instance)
        path = joinpath(current_dir, repo_path * instance)
        n,m,distances = readInstance_tsp(path)
        println("n : ", n)
        println("m : ", m)
        for p in test_files[instance]
            println("p : ", p)
            for k in k_list
                println("k : ", k)
                for method_init in methods_init
                    for method_mise_a_jour in methods_mise_a_jour
                        if method_mise_a_jour == "Aleatoire" && instance == "qa194.tsp"
                            continue
                        end
                        methods = [method_init, method_mise_a_jour]
                        println("methods : ", methods)
                        obj, time = algo(n,m,distances,p,k,methods)
                        write(out,instance*";"*string(n)*";"*string(m)*";"*string(p)*";"*string(k)*";"*string(method_init)*";"*string(method_mise_a_jour)*";"*string(obj)*";"*string(time)*"\n")                
                    end
                end
                flush(out)
            end
        end
    end
    close(out)
end


benchmark_grp8()

