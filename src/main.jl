import Pkg
Pkg.activate(".") # C'est mieux de pas polluer l'environnement de base

include("readData.jl")
include("groupes/groupe1.jl")
include("groupes/groupe2.jl")
include("groupes/groupe3.jl")
include("groupes/groupe4.jl")
include("groupes/groupe5.jl")
include("groupes/groupe6.jl")
include("groupes/groupe7.jl")
include("groupes/groupe8.jl")
include("groupes/groupe9.jl")
include("groupes/groupeA.jl")


function main_grp1(filepath)
    # Execution du travail du groupe 1 sur la formulation faible/forte et l'heuritique gloutonne
    println("Résolution pour le fichier : $filepath")
    n, m, opening_cost, cost_connection = read_data(filepath)
    main_pls(n, m, opening_cost, cost_connection)
    main_heurGlou(n, m, opening_cost, cost_connection)
end

function main_grp2(filepath)
    println("Résolution pour le fichier : $filepath")
    n, m, opening_cost, cost_connection = read_data(filepath)
    main_pls_bis(n, m, opening_cost, cost_connection)
end

function main_grp3(filepath)
    # Execution du travail du groupe 3 : Problème de p-centres
    println("Résolution pour le fichier : $filepath avec 5 centres")
    if endswith(filepath, ".tsp")
        n, m, distances = readInstance_tsp(filepath)
        main_p_centres(n, m, distances, 5)
    else
        n, m, _, cost_connection = read_data(filepath)
        main_p_centres(n, m, cost_connection, 5)
    end
end 

function main_grp4(filepath)
    # Execution du travail du groupe 4
    println("Résolution pour le fichier : $filepath")
    n, m, distances = read_data_npc(filepath)

    obj, temps, y, z = NPC(n,m,p,distances)
    println("Valeur distance max = ", obj, "\n \n")
    println("Sites ouverts = ", y, "\n \n")
    println("Temps de calcul = ", temps, "\n \n")
end

function main_grp5(filepath)
    # Execution du travail du groupe 5 : PLS Primal Dual
    println("Résolution pour le fichier : $filepath")
    n, m, opening_cost, cost_connection = read_data(filepath)
    main_pls_primal_dual(n, m, opening_cost, cost_connection)
end

function main_grp6(filepath)
    # Execution du travail du groupe 6
    println("Résolution pour le fichier : $filepath")
    n, m, opening_cost, cost_connection = read_data(filepath)
    p = 2
    println("main_stable: ",main_stable(n, m, cost_connection,p))
    println("resol_p_centre_set_cover: ",resol_p_centre_set_cover(cost_connection,p,maximum(cost_connection))[1])
end

function main_grp7(filepath)
    # Execution du travail du groupe 7
end

function main_grp8(filepath)
    # Execution du travail du groupe 8
end

function main_grp9(filepath)
    # Execution du travail du groupe 9 : Problème de p-centres connexe
    p = 5   
    F, R, distanceMatrix = readInstance_tsp(filepath)
    threshold = Int(floor(mean(distanceMatrix)))
    println("Résolution pour le fichier : $filepath avec $p centres et un seuil $threshold")
    println("Résolution de la formulation SCFF")
    SCFF(F, R, distanceMatrix, threshold, p, false)
    println("Résolution de la formulation MCFF")
    MCFF(F, R, distanceMatrix, threshold, p, false)
    println("Résolution de la formulation MTZ")
    MTZ(F, R, distanceMatrix, threshold, p, false)
end 

function main_grpA(filepath)
    # Execution du travail du groupe A
    println("Résolution pour le fichier : $filepath")
    n, m, distances = readInstance_tsp(filepath)
    p = 5
    algo_fixation_z(n, m, p, distances)
end

function main_grp10(filepath)
    # Execution du travail du groupe 10
    println("Résolution pour le fichier : $filepath")
    n, m, distance_matrix, Abs, Ord = readInstance_tsp_adapted(filepath)
    p = 2
    clusters = 16
    obj, y, time = p_center_clustering(n, m, distance_matrix, Abs, Ord, p, nb_clusters=clusters, npci=false, time_limit=180)
    println("Objective value: ", obj)
    println("Time: ", time)
    println("Open sites: ", findall(y -> y == 1, y))
end

function main()
    # Problème de localisation simple
    # # Exemple d'utilisation avec le fichier "Example.txt"
    # n, m, opening_cost, cost_connection = read_data("data/instTest.txt")

    # # Affichage des résultats
    # println("n = $n") # nombre de clients
    # println("m = $m") # nombres de sites
    # println("Opening Cost : $opening_cost")
    # println("Cost Connection : $cost_connection")
    # println()

    # # Execution du code du groupe 1
    # # main_grp1("data/instTest.txt")

    # Execution du code du groupe 2
    # main_grp2("data/instTest.txt")
    # benchmark_grp2()


    # Problème de P centres

    # n, m, distances = readInstance_tsp("MLA_PROJET/tsp_data/dj38.tsp")
    # println("n = m = $n") # nombre de sites
    # println("Distances : $opening_cost")

    # Execution du code du groupe 3
    # main_grp3("tsp_data/dj38.tsp")
    # benchmark_grp3()

    # # Execution du code du groupe 5
    # main_grp5("data/instTest.txt")

    # Execution du code du groupe 6
    # main_grp6("data/instRand_100_100_1.txt")

    # Execution du code du groupe 7

    # Execution du code du groupe 8

    # Execution du code du groupe 9

    # Execution du code du groupe A
    main_grpA("tsp_data/dj38.tsp")

    # Execution du code du groupe 10

end

main()
