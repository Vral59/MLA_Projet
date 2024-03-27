using JuMP
using CPLEX
using CSV
using DataFrames

# Mohamed Hamdi et Daria Stepanova
# formulation NPC





function read_data_npc(file_path)
    """Lis n,m et les distances"""
    
    file = open(file_path, "r")

    filename = readline(file)
    println("Lecture du fichier : $filename")

    n, m, _ = parse.(Int, split(readline(file)))

    cost_connection = zeros(Int, m, n)

    for i in 1:n
        line = split(readline(file))
        cost_connection[:, i] .= parse.(Int, line[3:end])
    end

    close(file)

    n, m = m,n

    return n, m, cost_connection
end



function distances_triées(n::Int, m::Int, distances::Matrix{Int})
    """Renvoit les distances différentes triées et leur nombre K"""
    # n nombre de clients
    # m nombre de sites
    
    distances_différentes = []
    for i in 1:n
        for j in 1:m
            if !(distances[i,j] in distances_différentes)
                push!(distances_différentes,distances[i,j])
            end
        end
    end
    
    distances_différentes = sort(distances_différentes)
    K = length(distances_différentes)

    return (K,distances_différentes)
end



function NPC(n::Int, m::Int, p::Int, distances::Matrix{Int}; silence::Bool = false)
    """formulation NPC"""
    # n nombre de clients
    # m nombre de sites
    # p nombre de sites à ouvrir
    # distances matrice des distances/coûts

    distances_K_D = distances_triées(n::Int, m::Int, distances)
    K, D = distances_K_D[1], distances_K_D[2]


    model = Model(CPLEX.Optimizer)
    if silence
        set_silent(model)
    end

    @variable(model, z[1:(K-1)], Bin)
    @variable(model, y[1:m], Bin)

    
    @constraint(model, [i in 1:n], sum(y[j] for j in 1:m)  == p)
    
    @constraint(model, [i in 1:n, k in 1:(K-1)], z[k]+sum(y[j] for j in 1:m if distances[i,j]<=D[k])  >= 1)



    @objective(model, Min, D[1]+sum((D[k]-D[k-1]) * z[k-1] for k in 2:K))

    optimize!(model)


    return objective_value(model), round(solve_time(model),digits=5), value.(y), value.(z)
end




function NPCi(n::Int, m::Int, p::Int, distances::Matrix{Int}; silence::Bool = false)
    """formulation NPC"""
    # n nombre de clients
    # m nombre de sites
    # p nombre de sites à ouvrir
    # distances matrice des distances/coûts

    distances_K_D = distances_triées(n::Int, m::Int, distances)
    K, D = distances_K_D[1], distances_K_D[2]


    model = Model(CPLEX.Optimizer)
    if silence
        set_silent(model)
    end

    @variable(model, z[1:K], Bin)
    @variable(model, y[1:m], Bin)

    
    @constraint(model, [i in 1:n], sum(y[j] for j in 1:m)  == p)
    
    @constraint(model, [i in 1:n, k in 1:K], z[k]+sum(y[j] for j in 1:m if distances[i,j]<=D[k])  >= 1)

    @constraint(model, [k in 1:(K-1)], z[k] >= z[k+1])



    @objective(model, Min, D[1]+sum((D[k]-D[k-1]) * z[k-1] for k in 2:K))

    optimize!(model)


    return objective_value(model), round(solve_time(model),digits=5), value.(y), value.(z)
end



function NPCir(n::Int, m::Int, p::Int, distances::Matrix{Int}; silence::Bool = false)
    """formulation NPC"""
    # n nombre de clients
    # m nombre de sites
    # p nombre de sites à ouvrir
    # distances matrice des distances/coûts

    distances_K_D = distances_triées(n::Int, m::Int, distances)
    K, D = distances_K_D[1], distances_K_D[2]


    model = Model(CPLEX.Optimizer)
    if silence
        set_silent(model)
    end

    @variable(model, z[1:K] >= 0)
    @variable(model, y[1:m], Bin)

    
    @constraint(model, [i in 1:n], sum(y[j] for j in 1:m)  == p)
    
    @constraint(model, [i in 1:n, k in 1:K], z[k]+sum(y[j] for j in 1:m if distances[i,j]<=D[k])  >= 1)

    @constraint(model, [k in 1:(K-1)], z[k] >= z[k+1])



    @objective(model, Min, D[1]+sum((D[k]-D[k-1]) * z[k-1] for k in 2:K))

    optimize!(model)


    return objective_value(model), round(solve_time(model),digits=5), value.(y), value.(z)
end




function main_npc(data, p)

    n, m, distances = read_data_npc(data)

    obj, temps, y, z = NPC(n,m,p,distances)
    println("Valeur distance max = ", obj, "\n \n")
    println("Sites ouverts = ", y, "\n \n")
    println("Temps de calcul = ", temps, "\n \n")

end



function lister_chemins_data(dossier::String="")
    """Renvoit la liste des chemins des fichiers de données"""
    # dossier: chemin du dossier 'data' 

    chemin_dossier = dossier == "" ? pwd() : joinpath(pwd(), dossier)
    
    println("Chemin du dossier : ", chemin_dossier) 
    
    chemins_data = String[]
    
    for fichier in readdir(chemin_dossier, join=true)
        println("Fichier trouvé : ", fichier) 
        if isfile(fichier)
            push!(chemins_data, fichier)
        end
    end
    
    return reverse(chemins_data)
end



function resultats_npc_csv(chemins_fichiers::Vector{String}, fichier_csv::String)
    """Sauvegarde les resultat dans un csv"""
    # chemins_fichiers: liste des chemins des fichiers
    # fichier_csv: chemin du fichier csv avec les resultats_csv

    df = CSV.read(fichier_csv, DataFrame)

    for chemin in chemins_fichiers
        nom_fichier = basename(chemin)
        n, m, distances = read_data_npc(chemin)
        obj, temps, y, z = NPC(n,m,3,distances)
        push!(df, (instances=nom_fichier, p=3, distance_max=obj, temps_calcul=temps))
    end

    CSV.write(fichier_csv, df)
end



function resultats_npci_csv(chemins_fichiers::Vector{String}, fichier_csv::String)
    """Sauvegarde les resultat dans un csv"""
    # chemins_fichiers: liste des chemins des fichiers
    # fichier_csv: chemin du fichier csv avec les resultats_csv

    df = CSV.read(fichier_csv, DataFrame)

    for chemin in chemins_fichiers
        nom_fichier = basename(chemin)
        n, m, distances = read_data_npc(chemin)
        obj, temps, y, z = NPCi(n,m,3,distances)
        push!(df, (instances=nom_fichier, p=3, distance_max=obj, temps_calcul=temps))
    end

    CSV.write(fichier_csv, df)
end



function resultats_npcir_csv(chemins_fichiers::Vector{String}, fichier_csv::String)
    """Sauvegarde les resultat dans un csv"""
    # chemins_fichiers: liste des chemins des fichiers
    # fichier_csv: chemin du fichier csv avec les resultats_csv

    df = CSV.read(fichier_csv, DataFrame)

    for chemin in chemins_fichiers
        nom_fichier = basename(chemin)
        n, m, distances = read_data_npc(chemin)
        obj, temps, y, z = NPCir(n,m,3,distances)
        push!(df, (instances=nom_fichier, p=3, distance_max=obj, temps_calcul=temps))
    end

    CSV.write(fichier_csv, df)
end




chemins_fichiers = lister_chemins_data("data")
fichier_csv = "resultats_npc_groupe4.csv" 
#resultats_npc_csv(chemins_fichiers, fichier_csv)
fichier_csv = "resultats_npci_groupe4.csv" 
#resultats_npci_csv(chemins_fichiers, fichier_csv)
fichier_csv = "resultats_npcir_groupe4.csv" 
#resultats_npcir_csv(chemins_fichiers, fichier_csv)



