# Mathis Azéma et Michel Sénégas
using JuMP
using CPLEX

function main_stable(n::Int, m::Int, cost_connection::Matrix{Int},p::Int)
    distances_K_D = distances_triées(n::Int, m::Int, cost_connection)
    K, D = distances_K_D[1], distances_K_D[2]
    kup=K
    klb=1
    while kup-klb>1
        delta= D[Int(floor((kup+klb)/2))]
        Arcs_Gp=creation_graphe_Gp(n,m ,cost_connection, delta)
        val=resolution_stable_max(n, Arcs_Gp)
        if val<=p
            kup=Int(floor((kup+klb)/2))
        else
            klb=Int(floor((kup+klb)/2))
        end
    end
    if kup==2
        delta= D[1]
        Arcs_Gp=creation_graphe_Gp(n,m ,cost_connection, delta)
        val=resolution_stable_max(n, Arcs_Gp)
        return D[1]
    else
        return D[kup]
    end
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

function resolution_stable_max(n::Int, Arcs_Gp::Vector{Any})
    """Calcul le stable dans le graphe Gp"""
    # n nombre de clients
    # Arcs G_p les arcs du graphe Gp
    
    model = Model(CPLEX.Optimizer)
    set_silent(model)

    @variable(model, x[1:n], Bin)

    
    @constraint(model, [(i, ip) in Arcs_Gp], x[i]+x[ip] <= 1)
    @objective(model, Max, sum(x[i] for i in 1:n))

    optimize!(model)

    return objective_value(model)
end

function creation_graphe_Gp(n::Int, m::Int, cost_connection::Matrix{Int}, delta::Int)
    """Calcul les arcs du graphe Gp."""
    # n nombre de clients, m le nombre de sites
    # delta distance pour créer le graphe G
    Arcs_G=[]
    sites_clients_delta=[[] for j in 1:m]
    for i in 1:n
        for j in 1:m
            if cost_connection[i,j]<= delta #Il existe un arc dans le graphe G si le site j est à moins de delta du client i
                push!(Arcs_G, (i,j))
                push!(sites_clients_delta[j], i)
            end
        end
    end
    Arcs_Gp=[]
    Arcs_Gp_exist=[[0 for ip in 1:n] for i in 1:n]
    for j in 1:m
        for i in sites_clients_delta[j]
            for ip in sites_clients_delta[j]
                if ip>i && Arcs_Gp_exist[i][ip]==0 # Il existe un arc dans Gp si (i,ip) ont un voisin commun dans le graphe G
                    push!(Arcs_Gp, (i,ip))
                    Arcs_Gp_exist[i][ip]=1
                end
            end
        end
    end
    return Arcs_Gp
end

"""Résout le problème de set cover.\\
   `relax::Bool` - Si l'on doit résoudre la relaxation continue uniquement"""
function set_cover(distances::Matrix{Int},δ::Int,relax=false)
    n,m = size(distances)

    model = Model(CPLEX.Optimizer)
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
    while LB + 1 < UB
        mean = (LB + UB) ÷ 2
        if f(mean)
            UB = mean
        else
            LB = mean
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
