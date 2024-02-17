"""
    read_data(file_path)

Lire les données à partir d'un fichier texte selon le format spécifié.

# Arguments
- `file_path::AbstractString`: Chemin d'accès au fichier texte.

# Returns
- `n::Int`: Nombre de lignes dans le fichier.
- `m::Int`: Nombre de valeurs dans chaque ligne.
- `opening_cost::Vector{Int}`: Vecteur contenant les deuxièmes chiffres de chaque ligne.
- `cost_connection::Matrix{Int}`: Matrice contenant les valeurs des lignes suivantes du fichier.

# Example
n, m, opening_cost, cost_connection = read_data("Example.txt")
"""
function read_data(file_path)
    # Ouvrir le fichier en mode lecture
    file = open(file_path, "r")

    # Lire la première ligne (nom du fichier)
    filename = readline(file)
    println("Lecture du fichier : $filename")

    # Lire la deuxième ligne pour obtenir n, m et ignorer le dernier élément (0)
    n, m, _ = parse.(Int, split(readline(file)))

    # Initialiser le vecteur et la matrice
    opening_cost = zeros(Int, n)
    cost_connection = zeros(Int, m, n)

    # Lire les lignes suivantes pour remplir le vecteur et la matrice
    for i in 1:n
        line = split(readline(file))
        opening_cost[i] = parse(Int, line[2])
        cost_connection[:, i] .= parse.(Int, line[3:end])
    end

    # Fermer le fichier
    close(file)

    n, m = m,n

    return n, m, opening_cost, cost_connection
end

