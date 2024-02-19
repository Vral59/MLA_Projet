include("readData.jl")
include("groupe1.jl")

function main_grp1(filepath)
    # Execution du travail du groupe 1 sur la formulation faible/forte et l'heuritique gloutonne
    println("Résolution pour le fichier : $filepath")
    n, m, opening_cost, cost_connection = read_data(filepath)
    main_pls(n, m, opening_cost, cost_connection)
    main_heurGlou(n, m, opening_cost, cost_connection)
end

function main_grp2(filepath)
    # Execution du travail du groupe 2
end 

function main_grp3(filepath)
    # Execution du travail du groupe 3
end 

function main()
    # Exemple d'utilisation avec le fichier "Example.txt"
    n, m, opening_cost, cost_connection = read_data("data/instTest.txt")

    # Affichage des résultats
    println("n = $n") # nombre de clients
    println("m = $m") # nombres de sites
    println("Opening Cost : $opening_cost")
    println("Cost Connection : $cost_connection")

    # Execution du code du groupe 1
    main_grp1("data/instTest.txt")

    # Execution du code du groupe 2

    # Execution du code du groupe 3
end

main()
