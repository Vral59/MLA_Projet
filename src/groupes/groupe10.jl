# Benjamin Jacquet et Christian El Maalouly

using JuMP
using CPLEX

# Fonction readInstance_tsp de readData, adaptée (retourne abs et ord en plus)
function readInstance_tsp_adapted(path::String)
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
                push!(Abs, x)
                push!(Ord, y)
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
    distance_matrix = zeros(Int64, n, n)
    for i in 1:n
        for j in 1:n
            distance_matrix[i, j] = round(Int64, sqrt((Abs[i] - Abs[j])^2 + (Ord[i] - Ord[j])^2))
        end
    end
    println("instance file $name parsed successfully")
    return n, m, distance_matrix, Abs, Ord
end

# Fonction du groupe 3 pour pc (adaptée pour le clustering)
function pc(I, m, cost_connection, nb_centres, time_limit=180)

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 8, "CPX_PARAM_TILIM" => time_limit))
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

    n = length(I)
    println("Resoudre le probleme du p-centre avec ", n, " clients et ", m, " sites ouverts, pour ", nb_centres, " centres.")

    @variable(model, x[i in I, 1:m] >= 0, Bin)
    @variable(model, y[1:m], Bin)
    @variable(model, r >= 0)

    @constraint(model, [i in I, j in 1:m], x[i, j] <= y[j])
    @constraint(model, [i in I], sum(x[i, :]) == 1)

    # p centres
    @constraint(model, sum(y) == nb_centres)

    # bound on r in PC1
    @constraint(model, [i in I, j in 1:m], cost_connection[i, j] * x[i, j] <= r)

    @objective(model, Min, r)

    optimize!(model)

    obj, x_val, y_val = 0, 0, 0
    if termination_status(model) == MOI.OPTIMAL
        obj = objective_value(model)
        x_val, y_val = value.(x), value.(y)
    end
    return obj, y_val, solve_time(model)
end

function distances_triées(n::Int, m::Int, distances::Matrix{Int})
    """Renvoit les distances différentes triées et leur nombre K"""
    # n nombre de clients
    # m nombre de sites

    distances_différentes = []
    for i in 1:n
        for j in 1:m
            if !(distances[i, j] in distances_différentes)
                push!(distances_différentes, distances[i, j])
            end
        end
    end

    distances_différentes = sort(distances_différentes)
    K = length(distances_différentes)

    return (K, distances_différentes)
end

function NPCi(I, m, cost_connection, nb_centres, time_limit=180)
    """formulation NPCi"""
    # I ensemble de clients
    # m nombre de sites
    # cost_connection matrice des distances/coûts
    # nb_centres nombre de sites à ouvrir
    n = length(I)
    new_costs = cost_connection[I, :]
    distances_K_D = distances_triées(n::Int, m::Int, new_costs)
    K, D = distances_K_D[1], distances_K_D[2]

    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 8, "CPX_PARAM_TILIM" => time_limit))
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output
    println("Resoudre le probleme du p-centre avec ", n, " clients et ", m, " sites ouverts, pour ", nb_centres, " centres.")

    function subset_K(i)
        subK = []
        for k in 1:K-1
            for j in 1:m
                if cost_connection[i, j] == D[k+1]
                    push!(subK, k)
                    break
                end
            end
        end
        return subK
    end

    @variable(model, z[1:K], Bin)
    @variable(model, y[1:m], Bin)

    @constraint(model, [i in I], sum(y[j] for j in 1:m) == nb_centres)
    @constraint(model, [i in I, k in subset_K(i)], z[k] + sum(y[j] for j in 1:m if cost_connection[i, j] <= D[k]) >= 1)
    # @constraint(model, [i in I, k in 1:K], z[k]+sum(y[j] for j in 1:m if cost_connection[i,j]<=D[k])  >= 1)
    @constraint(model, [k in 1:(K-1)], z[k] >= z[k+1])

    @objective(model, Min, D[1] + sum((D[k] - D[k-1]) * z[k-1] for k in 2:K))

    optimize!(model)

    obj, z_val, y_val = 0, 0, 0
    if termination_status(model) == MOI.OPTIMAL
        obj = objective_value(model)
        z_val, y_val = value.(z), value.(y)
    end
    return obj, y_val, solve_time(model)
end

# Fonction pour générer les clusters et leur representants initials
function generate_I(n, abs, ord, nb_clusters)

    grid_size = sqrt(nb_clusters)
    abs_range = (minimum(abs), maximum(abs))
    ord_range = (minimum(ord), maximum(ord))

    clusters = [[[] for _ in 1:grid_size] for _ in 1:grid_size]
    representatives = []

    abs_grid_range = (abs_range[2] - abs_range[1]) / grid_size
    ord_grid_range = (ord_range[2] - ord_range[1]) / grid_size

    for i in 1:n
        # Find the grid coordinates of the client
        if abs[i] == abs_range[2]
            grid_x = Int(grid_size)
        else
            grid_x = Int(floor((abs[i] - abs_range[1]) / abs_grid_range + 1))
        end
        if ord[i] == ord_range[2]
            grid_y = Int(grid_size)
        else
            grid_y = Int(floor((ord[i] - ord_range[1]) / ord_grid_range + 1))
        end

        # Add the coordinate to the corresponding cluster
        push!(clusters[grid_x][grid_y], i)
        if length(clusters[grid_x][grid_y]) == 1
            push!(representatives, i)
        end
    end

    return clusters, representatives
end

function p_center_clustering(n, m, cost_connection, abs, ord, nb_centres; npci=false, nb_clusters=100, time_limit=180)
    start_time = time()
    println("method: ", npci ? "NPCi" : "PC")

    LB = 0
    UB = Inf

    obj = 0
    best_y = zeros(Int, m)

    clusters, representatives = generate_I(n, abs, ord, nb_clusters)

    # Tant que la borne inférieure est inférieure à la borne supérieure et que le temps n'est pas écoulé
    while LB < UB && time() - start_time < time_limit
        # Résoudre le problem du p-center avec les représentants actuels
        if npci
            obj, y_val = NPCi(representatives, m, cost_connection, nb_centres, time_limit)
        else
            obj, y_val = pc(representatives, m, cost_connection, nb_centres, time_limit)
        end
        println("one iter pc done")
        open_sites = findall(y -> y == 1, y_val)
        println("obj: ", obj, " y_val: ", open_sites)

        # Si une solution est trouvée
        if obj != 0
            LB = obj
            value_for_all_clients = 0

            # Chercher dans chaque cluster le client le plus loin d'un site ouvert
            for row in clusters
                for cluster in row
                    farthest_distance = 0
                    farthest_client = 0

                    # Chercher parmit les clients du cluster sauf les representants
                    for i in cluster
                        if !(i in representatives)
                            # Chercher le site ouvert le plus proche du client
                            closest_site = cost_connection[i, open_sites[1]]
                            for j in open_sites
                                if cost_connection[i, j] < closest_site
                                    closest_site = cost_connection[i, j]
                                end
                            end

                            if closest_site > farthest_distance
                                farthest_distance = closest_site
                                farthest_client = i
                            end
                        end
                    end

                    # Ajouter le client le plus loin à la liste des representants (sauf si pas de client trouvé)
                    if farthest_client != 0
                        push!(representatives, farthest_client)

                        # Mettre à jour la valeur de la solution pour tous les clients
                        if farthest_distance > value_for_all_clients
                            value_for_all_clients = farthest_distance
                        end
                    end
                end
            end

            # Mettre à jour la borne supérieure si la solution est meilleure, ainsi que les sites de la meilleur solution
            if value_for_all_clients < UB
                if value_for_all_clients < LB
                    UB = LB
                else
                    UB = value_for_all_clients
                end
                best_y = copy(y_val)
            end

        else
            println("An iteration had no solution or took too long to solve")
            break
        end
        println("LB: ", LB, " UB: ", UB, " Time: ", time() - start_time)
    end

    return UB, best_y, time() - start_time


end
