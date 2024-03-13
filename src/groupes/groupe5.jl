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

    f_residual = [f[j] for j = 1:m]

    non_affected = n

    while non_affected > 0

        for i = 1:n

            if sum(x[i,:]) == 0 # non-affected client, check if condition 1. is satisfied
                for j = 1:m
                    if y[j] == 1 && v[i] == c[i,j]
                        x[i,j] = 1 
                    end
                end
                if sum(x[i,:]) >= 1 # if the client is now affected to one (or more) site
                    non_affected -=1
                end
            end

            if sum(x[i,:]) == 0  # still a non-affected client

                v[i] += 1  # Increase v[i]

                for j = 1:m 

                    # This block is useful only when f[j] = 0. Then,  3. is satisfied right away -> open j.
                    if f_residual[j] == 0
                        y[j] = 1
                        for i_prime = 1:n
                            if sum(x[i_prime,:]) == 0  && v[i_prime] >= c[i_prime,j]
                                x[i_prime, j] = 1
                                non_affected -=1
                                if non_affected == 0
                                    break       
                                end
                            end
                        end
                    end

                    if y[j] == 0 

                        if v[i] > c[i,j] && f_residual[j] >= 1 # Condition 2. implies than w[i,j] should be incremented

                            w[i,j] += 1
                            f_residual[j] -= 1

                            if f_residual[j] == 0 # Condition 3. is satisfied, open j
                                y[j] = 1
                                for i_prime = 1:n
                                    if sum(x[i_prime,:]) == 0  && v[i_prime] >= c[i_prime,j] 
                                        x[i_prime, j] = 1
                                        non_affected -=1
                                        if non_affected == 0
                                            break       
                                        end
                                    end
                                end
                            end 

                        end 
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