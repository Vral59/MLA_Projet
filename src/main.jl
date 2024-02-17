include("readData.jl")


function main()
    # Exemple d'utilisation avec le fichier "Example.txt"
    n, m, opening_cost, cost_connection = read_data("data/instTest.txt")

    # Affichage des résultats
    println("n = $n") # nombre de clients
    println("m = $m") # nombres de sites
    println("Opening Cost : $opening_cost")
    println("Cost Connection : $cost_connection")
end


main()