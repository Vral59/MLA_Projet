# Ines Chergui et Margot Boyer
using JuMP
using CPLEX
using Dates
using Statistics
using StatsBase

include("../readData.jl")

current_dir = pwd()

path = joinpath(current_dir, "MLA_Projet\\tsp_data\\dj38.tsp")
#path = "MLA_Projet\\data\\ga250a-4.txt"


n,m,distances = readInstance_tsp(path)
println("n : ", n)
println("m : ", m)

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
    @constraint(model, [i in I], sum(x[i,:])  == 1)
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

    #n_variables = num_variables(model)
    #n_constraints = num_constraints(model; count_variable_in_set_constraints=false) # Count only structural constraints

    status = termination_status(model)
    @assert status == MOI.OPTIMAL || status == MOI.TIME_LIMIT

    obj, x_val = missing, missing
    if primal_status(model) == MOI.FEASIBLE_POINT
        obj = objective_value(model)
        x_val = value.(x)
    end
    return obj, x_val
end

function mise_a_jour(n::Int, m::Int, I::Vector{Int}, x_val::Matrix{Float64}, distances::Matrix{Int}, method::String = "Distance_max")
    """
    Met à jour le sous-ensemble de clients selon différentes méthodes
    """
    if method == "Distance_max"
        """
        Ajoute le client qui maximise la distance aux usines ouvertes courantes
        """
        r = x_val .* distances
        r[I,:] .= -1
        client = argmax(r)[1]

    elseif method == "Aleatoire"
        """
        Ajoute un client aléatoirement
        """
        Ibar = collect(setdiff(1:n, I))
        client = rand(Ibar)
    end
    push!(I,client)
    return I
end


function initialise_I(n::Int,m::Int,distances::Matrix{Int},k::Int=1,method="Aleatoire")
    """ Initialise le sous-ensemble de clients I
    Args : 
        - k : nombre de clients de l'ensemble d'inialisation I
    """
    if method == "Distance_moyenne"
        distances_moyenne = mean(distances, dims =  1)
        println("distance moyenne : ", distances_moyenne)
        
        clients_tries = sortperm(vec(distances_moyenne), rev=true)
        I = Vector{Int}(clients_tries[1:k])
    elseif method == "Aleatoire"
        #I = Vector{Int}([rand(1:n)])
        I = sample(1:n,k,replace = false)
    end
    return I
end

function algo(n::Int, m::Int, distances::Matrix{Int}, p::Int)
    #I = Vector{Int}([rand(1:n)])
    start_time = now()
    I = initialise_I(n,m,distances,2, "Distance_moyenne")
    println("I : ", I)
    lb = 0
    ub = Inf
    while ub > lb
        println("Difference ub et lb : ", ub - lb)
        vopt, time, node_count, x_val, y_val, n_variables, n_constraints, lb = pc(p,m,I,distances)
        lb = vopt
        println("Temps passe solve pc : ",time )
        ry, x_val = rayon(n,m,distances,y_val)
        ub = min(ub, ry)
        I = mise_a_jour(n, m, I, x_val, distances, "Aleatoire")
        println("Nouvel ensemble I : ",I)
        println("Nouvelle borne inf : ", lb)
        println("Nouvelle borne sup : ", ub)
    end
    execution_time = Dates.value(now() - start_time)/1000
    println("Temps d'exécution(s) : ", execution_time)
end


algo(n,m,distances,5)

