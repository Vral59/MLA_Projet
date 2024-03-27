# Mathis Azéma et Michel Sénégas
using JuMP
using CPLEX

function main_stable(n::Int, m::Int, cost_connection::Matrix{Int},p::Int)
    """  Recherche dichotomique avec la méthode du stable """
    distances_K_D = distances_triées(n::Int, m::Int, cost_connection)
    K, D = distances_K_D[1], distances_K_D[2]
    kup=K
    klb=1
    solutions=[[0.0] for d in D] #Stocke les solutions du stable maximum pour une valeur de d -> utile pour la construction d'une solution réalisable
    while kup-klb>1 #Actualisation des bornes de la dichotomie
        ind=Int(floor((kup+klb)/2))
        delta= D[ind]
        Arcs_Gp=creation_graphe_Gp(n,m ,cost_connection, delta)
        val, sol_delta=resolution_stable_max(n, Arcs_Gp)
        if val<=-1 #val = -1 signifie que le stable maximum n'a pas été trouvé dans le temps limite
            return (-1, -1)
        end
        solutions[ind]=sol_delta
        if val<=p
            kup=Int(floor((kup+klb)/2))
        else
            klb=Int(floor((kup+klb)/2))
        end
    end
    if kup==2 #Gestion du cas où kup==2
        delta= D[1]
        Arcs_Gp=creation_graphe_Gp(n,m ,cost_connection, delta)
        val, sol=resolution_stable_max(n, Arcs_Gp)
        if val<=-1
            return (-1, -1)
        end
        solutions[1]=sol
        if val<=p
            val_sol_realisable=build_solution_after_stable(n, m, p, cost_connection, D[1], solutions[1])
            return D[1], val_sol_realisable
        else
            if solutions[2]==[0.0]
                Arcs_Gp=creation_graphe_Gp(n,m ,cost_connection, D[2])
                val, sol=resolution_stable_max(n, Arcs_Gp)
                if val<=-1
                    return (-1, -1)
                end
                solutions[2]=sol
            end
            Arcs_Gp=creation_graphe_Gp(n,m ,cost_connection, D[2])
            val_sol_realisable=build_solution_after_stable(n, m, p, cost_connection, D[kup], solutions[kup])
            return D[2], val_sol_realisable
        end
    else
        if solutions[kup]==[0.0]
            Arcs_Gp=creation_graphe_Gp(n,m ,cost_connection, D[kup])
            val, sol=resolution_stable_max(n, Arcs_Gp)
            if val<=-1
                return (-1, -1)
            end
            solutions[kup]=sol
        end
        val_sol_realisable=build_solution_after_stable(n, m, p, cost_connection, D[kup], solutions[kup])
        return D[kup], val_sol_realisable
    end
end 

function build_solution_after_stable(n::Int, m::Int, p::Int, cost_connection::Matrix{Int}, delta::Int, solution_stable::Vector{Float64})
    """Calcul d'une solution réalisable à partir d'un stable maximum
        On ouvre d'abord des sites liés au client du stable puis on ajoute les clients non affectés à leur site ouvert le plus proche"""
    client_stable=[i for i in 1:n if solution_stable[i]>0.1]
    client_affected=[0 for i in 1:n]
    site_is_opened=[0 for j in 1:m]
    solution=[0 for i in 1:n]
    nb_sites_ouverts=0
    while nb_sites_ouverts < length(client_stable)
        current_client=nothing
        for id in client_stable
            if client_affected[id]==0
                current_client=id
                break
            end
        end
        possible_site=[j for j in 1:m if site_is_opened[j]==0]
        open_site=possible_site[argmin([cost_connection[current_client, j] for j in 1:m if site_is_opened[j]==0])] #On ouvre le site le plus proche au client courant présent dans le stable
        nb_sites_ouverts+=1
        site_is_opened[open_site]=1
        client_affected[current_client]=1
        solution[current_client]=open_site
        for i in 1:n
            if client_affected[i]==0 && cost_connection[i, open_site]<=delta
                solution[i]=open_site
                client_affected[i]=1
            end
        end
    end
    for i in 1:n #Affectation des clients restants
        if client_affected[i]==0
            client_affected[i]=1
            possible_site=[j for j in 1:m if site_is_opened[j]==1]
            solution[i]=possible_site[argmin([cost_connection[i,j] for j in 1:m if site_is_opened[j]==1])]
        end
    end
    return maximum([cost_connection[i, solution[i]] for i in 1:n])
end

function distances_triées(n::Int, m::Int, distances::Matrix{Int})
    """Renvoit les distances différentes triées et leur nombre K"""
    # n nombre de clients
    # m nombre de sites
    
    distances_differentes_ = Set(distances)

    distances_differentes = [d for d in distances_differentes_]
    
    distances_differentes = sort(distances_differentes)
    K = length(distances_differentes)

    return (K,distances_differentes)
end

function resolution_stable_max(n::Int, Arcs_Gp::Vector{Tuple{Int64, Int64}})
    """Calcul le stable dans le graphe Gp"""
    # n nombre de clients
    # Arcs G_p les arcs du graphe Gp
    
    model = Model(CPLEX.Optimizer)
    set_silent(model)
    set_time_limit_sec(model, 30)

    @variable(model, x[1:n], Bin)

    
    @constraint(model, [(i, ip) in Arcs_Gp], x[i]+x[ip] <= 1)
    @objective(model, Max, sum(x[i] for i in 1:n))

    optimize!(model)

    if primal_status(model) == MOI.FEASIBLE_POINT
        if termination_status(model) == MOI.TIME_LIMIT
            return (-1, nothing)
        else 
            return objective_value(model), JuMP.value.(x)
        end
    else
        return (-1, nothing)
    end
end

function creation_graphe_Gp(n::Int, m::Int, cost_connection::Matrix{Int}, delta::Int)
    """Calcul les arcs du graphe Gp."""
    binary_matrix=cost_connection.<=delta
    link_matrix= binary_matrix*transpose(binary_matrix)
    for i in 1:n
        link_matrix[i,i]=0
    end
    index=findall(x -> x > 0, link_matrix)
    Arcs_Gp=[Tuple(arc) for arc in index]
    return Arcs_Gp
end

"""Résout le problème de set cover.\\
   `relax::Bool` - Si l'on doit résoudre la relaxation continue uniquement"""
function set_cover(distances::Matrix{Int},δ::Int,relax=false)
    n,m = size(distances)

    model = Model(CPLEX.Optimizer)
    set_time_limit_sec(model,10)
    set_silent(model)

    @variable(model,y[1:m],Bin)

    @constraint(model,[i in 1:n],sum(y[j] for j in 1:m if distances[i,j] ≤ δ) ≥ 1)

    @objective(model,MIN_SENSE,sum(y))

    if relax
        relax_integrality(model)
    end

    optimize!(model)

    if primal_status(model) == FEASIBLE_POINT
        return objective_value(model),value.(y)
    else
        return Inf,zeros(m)
    end
end

"""Heuristique d'arrondi simple pour le problème de set cover."""
function set_cover_arrondi(distances::Matrix{Int},δ::Int)
    v,y = set_cover(distances,δ,true)
    if isinf(v)
        return v,y
    end

    n,m = size(distances)
    argmaxy = ones(Int,n)
    for i in 1:n
        for j in 1:m
            if distances[i,j] ≤ δ
                if y[j] > y[argmaxy[i]]
                    argmaxy[i] = j
                end
            end
        end
        y[argmaxy[i]] = 1
    end
    replace!(s -> s < 1-2^-35 ? 0 : s,y)
    return sum(y),y
end

"""Heuristique d'arrondi complexe pour le problème de set cover."""
function set_cover_arrondi_mixte(distances::Matrix{Int},δ::Int)
    _,y = set_cover(distances,δ,true)

    n,m = size(distances)
    A = falses(n,m)
    setsize = zeros(Int,n)
    for j in 1:m
        for i in 1:n
            if distances[i,j] ≤ δ
                A[i,j] = true
                setsize[i] += 1
            end
        end
    end
    b = ones(n)

    i = argmin(ii -> setsize[ii] / b[ii],1:n)
    if setsize[i] == 0
        return Inf,y
    end
    while b[i] > 0
        j_order = sortperm(y,rev=true)
        for j in j_order
            if A[i,j]
                y[j] = b[i]
                v = b[i]
                for ii in 1:n
                    if A[ii,j]
                        A[ii,j] = false
                        setsize[ii] -= 1
                        b[ii] -= v
                        if b[ii] < 0 b[ii] = 0 end
                    end
                end
            end
        end
        replace!(s -> s == 0 ? 1 : s,setsize)
        i = argmin(ii -> setsize[ii] / b[ii],1:n)
    end
    replace!(s -> s < 1-2^-35 ? 0 : s,y)
    return sum(y),y
end

"""Heuristique pour le problème de set cover en O(mn²)"""
function set_cover_heuristique(distances::Matrix{Int},δ::Int)
    n,m = size(distances)
    A = falses(n,m)
    setsize = zeros(Int,n)
    for j in 1:m
        for i in 1:n
            if distances[i,j] ≤ δ
                A[i,j] = true
                setsize[i] += 1
            end
        end
    end
    b = ones(n)
    y = zeros(m)

    i = argmin(ii -> setsize[ii] / b[ii],1:n)
    if setsize[i] == 0
        return Inf,y
    end
    while b[i] > 0
        j_order = sortperm(y)
        for j in j_order
            if A[i,j]
                y[j] = b[i]
                v = b[i]
                for ii in 1:n
                    if A[ii,j]
                        A[ii,j] = false
                        setsize[ii] -= 1
                        b[ii] -= v
                        if b[ii] < 0 b[ii] = 0 end
                    end
                end
            end
        end
        replace!(s -> s == 0 ? 1 : s,setsize)
        i = argmin(ii -> setsize[ii] / b[ii],1:n)
    end
    return sum(y),y
end

"""Calcule la plus petite valeur telle que la fonction en entrée renvoie `true`\\
   `f : Int -> Bool` - Fonction booléenne **croissante** à évaluer."""
function dichotomie(f::Function,LB::Int,UB::Int)
    t0 = time()
    while LB + 1 < UB
        mean = (LB + UB) ÷ 2
        if f(mean)
            UB = mean
        else
            LB = mean
        end
        if time()-t0>30
            return Int(1e12)
        end
    end
    return UB
end

"""Résout le problème de p-centre à l'aide d'une série de problèmes de set cover."""
function resol_p_centre_set_cover(distances::Matrix{Int},p::Int,UB::Int)
    # Phase 1
    UB = dichotomie(δ -> set_cover_arrondi(distances,δ)[1] ≤ p,0,UB)
    LB = dichotomie(δ -> set_cover(distances,δ,true)[1] ≤ p,0,UB)
    # Phase 2
    OPT = dichotomie(δ -> set_cover(distances,δ)[1] ≤ p,LB,UB)
    y = set_cover(distances,OPT)[2]
    return OPT,y
end

############ REMOVE BELOW

function resol_p_centre_set_cover_bounds(distances::Matrix{Int},p::Int,UB::Int,phase2::Bool)
    # Phase 1
    UB = dichotomie(δ -> set_cover_arrondi(distances,δ)[1] ≤ p,0,UB)
    LB = dichotomie(δ -> set_cover(distances,δ,true)[1] ≤ p,0,UB)
    if !phase2
        return LB,UB,missing
    end
    # Phase 2
    OPT = dichotomie(δ -> set_cover(distances,δ)[1] ≤ p,LB,UB)
    return LB,UB,OPT
end
