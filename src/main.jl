include("readData.jl")
include("groupes/groupe1.jl")
include("groupes/groupe2.jl")
include("groupes/groupe3.jl")


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
    println("Résolution pour le fichier : $filepath")
    n, m, opening_cost, cost_connection = read_data(filepath)
    main_p_centres(n, m, cost_connection, 3)
end 

function main_grp4(filepath)
    # Execution du travail du groupe 4
end

function main_grp5(filepath)
    # Execution du travail du groupe 5
end

function main_grp6(filepath)
    # Execution du travail du groupe 6
end


function main()
    # Exemple d'utilisation avec le fichier "Example.txt"
    n, m, opening_cost, cost_connection = read_data("data/instTest.txt")

    # Affichage des résultats
    println("n = $n") # nombre de clients
    println("m = $m") # nombres de sites
    println("Opening Cost : $opening_cost")
    println("Cost Connection : $cost_connection")
    println()

    # Execution du code du groupe 1
    # main_grp1("data/instTest.txt")

    # Execution du code du groupe 2
    # main_grp2("data/instTest.txt")
    benchmark_grp2()

    # Execution du code du groupe 3
    main_grp3("data/instRand_50_50_1.txt")
end

main()
