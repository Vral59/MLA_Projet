#Arthur Divanovic et Axel Navarro

function main_pls_primal_dual(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int})
    obj, x, y, resolution_time = PLS_primal_dual(n, m, opening_cost, cost_connection)

    println("Coût PLS Primal-Dual= ", obj)

    println("Solutions (site j : [clients i affectés à j])")
    display_solution(n,x,y)

    println("Temps de résolution : ", resolution_time)
end 

function PLS_primal_dual(n::Int, m::Int, opening_cost::Vector{Int}, cost_connection::Matrix{Int})
    f = opening_cost
    c = cost_connection

    start_time = time()
    v, w, x, y = phase1(n,m,f,c)
    x, y = phase2(n,c,w,x,y)
    obj = sum(x .* c) + sum(y .* f)
    end_time = time() - start_time
    return obj, x, y, end_time
end

function phase1(n::Int, m::Int, f::Vector{Int}, c::Matrix{Int})::Tuple{Vector{Int},Matrix{Int},Matrix{Int},Vector{Int}}
    v = zeros(Int, n)
    w = zeros(Int, n, m)
    x = zeros(Int, n, m)
    y = zeros(Int, m)

    # Matrix to identify the indices i,j s.t. v_i = w_i,j should be maintained.
    relation_maintained = zeros(Int, n, m)

    all_affected = false
    while !all_affected

        # 1. v_i = c_i,j and facility j is open.
        for i = 1:n
            for j = 1:m
                if v[i] == c[i,j] && y[j] == 1
                    # Affect i to j.
                    x[i,j] = 1
                end
            end
        end

        # 2. v_i = c_i,j and facility j is closed.
        for i = 1:n
            for j = 1:m
                if v[i] == c[i,j] && y[j] == 0
                    # Maintain v_i = w_i,j.
                    relation_maintained[i,j] = 1
                end
            end
        end

        # 3. sum(w_i,j) = f_j for facility j.
        for j = 1:m
            if sum(w[:,j]) == f[j]
                # Open facility j.
                y[j] = 1

                for i = 1:n
                    if sum(x[i,:]) == 0 && v[i] >= c[i,j]
                        # Affect the non-affected clients s.t. v_i >= c_i,j.
                        x[i,j] = 1
                    end
                end

            end
        end

        # Check if all clients are affected and modify v if not.
        all_affected = true
        for i = 1:n
            if sum(x[i,:]) == 0
                # Client i is not affected.
                all_affected = false
                # Increase v_i.
                v[i] += 1
                for j = 1:m 
                    # Increase w_i,j for all the facilities j considered in step 2.
                    if relation_maintained[i,j] == 1
                        w[i,j] += 1
                    end
                end
            end
        end
    end

    return v, w, x, y

end

function phase2(n::Int, c::Matrix{Int}, w::Matrix{Int}, x::Matrix{Int}, y::Vector{Int})
    V = findall(x -> x == 1, y)
    adj = zeros(Int, length(V), length(V))

    #Construction de la matrice d'adjacence
    for idx1 = 1:length(V)
        for idx2 = idx1+1:length(V)
            j1 = V[idx1]
            j2 = V[idx2]
            for i = 1:n
                if w[i,j1] > 0 && w[i,j2] > 0
                    adj[idx1,idx2] = 1
                    adj[idx2, idx1] = 1
                end
            end
        end
    end

    #Index des sommets du stable X
    X_idx = zeros(Int, length(V))

    while true
        # Trouver le sommet de V \ X avec le moins de voisins dans X
        best_idx = argmin([X_idx[idx] == 1 ? typemax(Int) : sum(adj[idx,:] .* X_idx) for idx in 1:length(V)])
        
        # Si ce sommet est relié à X ou si tous les sommets ont été ajouté -> arrêter
        if sum(adj[best_idx, :] .* X_idx) > 0 || sum(X_idx) == length(V)
            break
        end
        
        # Sinon jouter le sommet à X
        X_idx[best_idx] = 1
    end

    V_menus_X_idx = [1-v for v in X_idx]
    V_menus_X = V[findall(x -> x == 1, V_menus_X_idx)]
    X = V[findall(x -> x == 1, X_idx)]

    
    for idx = 1:length(V)
        j = V[idx]
    
        # Si j appartient à V \ X
        if V_menus_X_idx[idx] == 1
            
    
            #Fermeture du site j
            y[j] = 0
    
            #Pour chaque client i 
            for i = 1:n
                
                # Si le client i était affecté à j
                if x[i,j] == 1
                    
                    # i n'est plus affecté à j 
                    x[i,j] = 0
                    
                    #Si i n'a plus de site affecté, on cherche un site j' dans X
                    if sum(x[i,:]) == 0
    
                        # On cherche le meilleur site j' appartenant X et voisin de j dans G
                        best_new_j = nothing
                        best_cost = typemax(Int)
        
                        for new_idx = 1:length(V)
                            # Si le sommet appartient à X, est voisin de j dans le graphe et a un meilleur coût 
                            if X_idx[new_idx] == 1 && adj[idx,new_idx] == 1 && c[i,V[new_idx]] < best_cost
                                best_new_j = V[new_idx]
                                best_cost = c[i,best_new_j]
                            end
                        end

                        # i est affecté à j'
                        x[i,best_new_j] = 1
                    end
                end
            end
        end
    end

    return x, y
end

function display_solution(n::Int, x::Matrix{Int}, y::Vector{Int})
    opens = findall(x -> x == 1, y)

    for j in opens
        clients = Int[]
        for i = 1:n
            if x[i,j] == 1
                push!(clients, i)
            end
        end
        println("site ", j, " : ", clients)
    end
end