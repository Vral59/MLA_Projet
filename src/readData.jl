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

function readInstance_tsp(path::String)
    """
    The parsing function for the TSP instance file
    Encoding options are not read
    """
    name = ""
    wtype = "round"
    nnodes = 0
    Abs = []
    Ord = []
    start_coords = false
    open(path) do f
        lines = readlines(f)
        for line in lines
            sline = split(line)
            if sline[1] == "EOF"
                break
            elseif start_coords
                p = parse(Int64, sline[1])
                x = parse(Float64, sline[2])
                y = parse(Float64, sline[3])
                Abs=vcat(Abs,x)
                Ord=vcat(Ord,y)
            elseif sline[1] == "NAME"
                name = sline[3]
            elseif sline[1] == "DIMENSION"
                nnodes = parse(Int64, sline[3])
            elseif sline[1] == "EDGE_WEIGHT_TYPE"
                if sline[3] == "EUC_2D"
                    wtype = "round"
                elseif sline[3] == "CEIL_2D"
                    wtype = "ceil"
                elseif sline[3] == "GEOM"
                    wtype = "geom"
                end
            elseif sline[1] == "NODE_COORD_SECTION"
                start_coords = true
            end
        end
    end
    @assert length(Abs) == length(Ord)
    n = length(Abs)
    m = n
    distance_matrix = zeros(Int64,n,n)
    for i in 1:n
        for j in 1:n
            distance_matrix[i,j] = round(Int64,sqrt((Abs[i]-Abs[j])^2+(Ord[i]-Ord[j])^2))
        end
    end
    println("instance file $name parsed successfully")
    return n, m, distance_matrix
end


