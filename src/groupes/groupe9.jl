# Salma Hammani et Francesco Gallo
# Problème du p-centre connexe
using JuMP
using CPLEX
using Statistics

include("../readData.jl")

function SCFF(F::Int, R::Int, distanceMatrix::Matrix{Int}, threshold::Int, p::Int, verbose=true, connected=true, relax=true)
    m = Model(CPLEX.Optimizer)
    if !verbose
        set_silent(m)
    end
    @variable(m, 0 <= xF[1:F, 1:F] <= 1) # x^F_ij avec i dans F et j and F
    @variable(m, 0 <= xr[1:F] <= 1) # x^F_ri avec i dans F et r la racine
    @variable(m, 0 <= y[1:F] <= 1) # y_i avec i dans F
    @variable(m, yr, Bin) # y_r où r est la racine
    @variable(m, 0 <= xR[1:F, 1:R] <= 1) # x^R_ij avec i dans F et j dans R
    @variable(m, r >= 0) # le rayon maximal
    @variable(m, f[1:F, 1:F] >= 0) # f_ij avec i et j dans F
    @variable(m, fr[1:F] >= 0) # f_ri avec i dans F

    @constraint(m, [j in 1:R], sum(xR[i, j] for i in 1:F) == 1)
    @constraint(m, [j in 1:R, i in 1:F], xR[i, j] <= y[i])
    @constraint(m, sum(y[i] for i in 1:F) == p)
    @constraint(m, [j in 1:R], r >= sum(distanceMatrix[i, j] * xR[i, j] for i in 1:F))

    @constraint(m, yr == 1)
    if connected
        for i in 1:F
            for j in 1:F
                if distanceMatrix[i, j] > threshold
                    @constraint(m, xF[i, j] == 0)
                end
            end
        end
        @constraint(m, [i in 1:F, j in 1:F], xF[i, j] <= y[i])
        @constraint(m, [i in 1:F, j in 1:F], xF[i, j] <= y[j])
        @constraint(m, [j in 1:F], sum(f[i, j] for i in 1:F) + fr[j] - sum(f[j, i] for i in 1:F) == y[j])
        @constraint(m, sum(fr[i] for i in 1:F) == sum(y[i] for i in 1:F))
        @constraint(m, [i in 1:F, j in 1:F], f[i, j] <= F * xF[i, j])
        @constraint(m, [i in 1:F], fr[i] <= F * xr[i])
        @constraint(m, sum(xr[i] for i in 1:F) == 1)
    end

    @objective(m, Min, r)

    if relax
        optimize!(m)
        relaxation = JuMP.objective_value(m)
    end

    for i in 1:F
        for j in 1:F
            set_binary(xF[i, j])
        end
        for j in 1:R
            set_binary(xR[i, j])
        end
        set_binary(y[i])
        set_binary(xr[i])
    end

    set_attribute(m, "CPXPARAM_TimeLimit", 15 * 60)
    start_time = time()
    optimize!(m)
    end_time = time()

    if !relax
        relaxation = JuMP.objective_bound(m)
    end

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        return relaxation, JuMP.objective_value(m), JuMP.objective_bound(m), JuMP.node_count(m), end_time - start_time
    else
        return relaxation, false, JuMP.objective_bound(m), JuMP.node_count(m), end_time - start_time
    end
end

function MCFF(F::Int, R::Int, distanceMatrix::Matrix{Int}, threshold::Int, p::Int, verbose=true, connected=true, relax=true)
    m = Model(CPLEX.Optimizer)
    if !verbose
        set_silent(m)
    end
    @variable(m, 0 <= xF[1:F, 1:F] <= 1) # x^F_ij avec i dans F et j and F
    @variable(m, 0 <= xr[1:F] <= 1) # x^F_ri avec i dans F et r la racine
    @variable(m, 0 <= y[1:F] <= 1) # y_i avec i dans F
    @variable(m, yr, Bin) # y_r où r est la racine
    @variable(m, 0 <= xR[1:F, 1:R] <= 1) # x^R_ij avec i dans F et j dans R
    @variable(m, r >= 0) # le rayon maximal
    @variable(m, fk[1:F, 1:F,1:F] >= 0) # f_k_ij avec i et j dans F
    @variable(m, frk[1:F,1:F] >= 0) # fk_ri avec i dans F


    @constraint(m, [j in 1:R], sum(xR[i, j] for i in 1:F) == 1)
    @constraint(m, [j in 1:R, i in 1:F], xR[i, j] <= y[i])
    @constraint(m, sum(y[i] for i in 1:F) == p)
    @constraint(m, [j in 1:R], r >= sum(distanceMatrix[i, j] * xR[i, j] for i in 1:F))

    @constraint(m, yr == 1)

    if connected
        for i in 1:F
            for j in 1:F
                if distanceMatrix[i, j] > threshold
                    @constraint(m, xF[i, j] == 0)
                end
            end
        end
        @constraint(m, [i in 1:F, j in 1:F], xF[i, j] <= y[i])
        @constraint(m, [i in 1:F, j in 1:F], xF[i, j] <= y[j])
        @constraint(m, [k in 1:F, j in 1:F ; k != j] , sum(fk[k, i, j] for i in 1:F) + frk[k,j] - sum(fk[k, j, i] for i in 1:F) == 0)
        @constraint(m, [k in 1:F] , sum(fk[k, i, k] for i in 1:F) + frk[k,k] - sum(fk[k, k, i] for i in 1:F) == y[k])
        @constraint(m, [k in 1:F] ,sum(frk[k,i] for i in 1:F) == y[k]) 
        @constraint(m, [k in 1:F ,i in 1:F, j in 1:F], fk[k, i, j] <= xF[i, j])
        @constraint(m, [k in 1:F ,i in 1:F], frk[k, i] <= xr[i])
        @constraint(m, sum(xr[i] for i in 1:F) == 1)
    end

    @objective(m, Min, r)

    if relax
        optimize!(m)
        relaxation = JuMP.objective_value(m)
    end

    for i in 1:F
        for j in 1:F
            set_binary(xF[i, j])
        end
        for j in 1:R
            set_binary(xR[i, j])
        end
        set_binary(y[i])
        set_binary(xr[i])
    end

    set_attribute(m, "CPXPARAM_TimeLimit", 15 * 60)
    start_time = time()
    optimize!(m)
    end_time = time()

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if !relax
        relaxation = JuMP.objective_bound(m)
    end

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        return relaxation, JuMP.objective_value(m), JuMP.objective_bound(m), JuMP.node_count(m), end_time - start_time
    else
        return relaxation, false, JuMP.objective_bound(m), JuMP.node_count(m), end_time - start_time
    end
end


function MTZ(F::Int, R::Int, distanceMatrix::Matrix{Int}, threshold::Int, p::Int, verbose=true, connected=true, relax=true)
    m = Model(CPLEX.Optimizer)
    if !verbose
        set_silent(m)
    end
    @variable(m, 0 <= xF[1:F, 1:F] <= 1) # x^F_ij avec i dans F et j and F
    @variable(m, 0 <= xr[1:F] <= 1) # x^F_ri avec i dans F et r la racine
    @variable(m, 0 <= y[1:F] <= 1) # y_i avec i dans F
    @variable(m, yr, Bin) # y_r où r est la racine
    @variable(m, 0 <= xR[1:F, 1:R] <= 1) # x^R_ij avec i dans F et j dans R
    @variable(m, r >= 0) # le rayon maximal
    @variable(m, u[1:F] >= 0) 
    @variable(m, ur >= 0) 

    @constraint(m, [j in 1:R], sum(xR[i, j] for i in 1:F) == 1)
    @constraint(m, [j in 1:R, i in 1:F], xR[i, j] <= y[i])
    @constraint(m, sum(y[i] for i in 1:F) == p)
    @constraint(m, [j in 1:R], r >= sum(distanceMatrix[i, j] * xR[i, j] for i in 1:F))

    @constraint(m, yr == 1)
    if connected
        for i in 1:F
            for j in 1:F
                if distanceMatrix[i, j] > threshold
                    @constraint(m, xF[i, j] == 0)
                end
            end
        end
        @constraint(m, [i in 1:F, j in 1:F], xF[i, j] <= y[i])
        @constraint(m, [i in 1:F, j in 1:F], xF[i, j] <= y[j])
        @constraint(m, [j in 1:F, k in 1:R], sum(xF[i, j] for i in 1:F) + xr[j] >= xR[j,k])
        @constraint(m, [k in 1:F, j in 1:F, j!= k], sum(xF[i, j] for i in 1:F) + xr[j] >= xF[j, k])
        @constraint(m, [i in 1:F,j in 1:F], (F+1) * xF[i,j] + u[i] <= u[j] + F )
        @constraint(m, [j in 1:F], (F+1) * xr[j] + ur <= u[j] + F )
        @constraint(m, sum(xr[i] for i in 1:F) == 1)
        @constraint(m, ur == 0)
    end

    @objective(m, Min, r)

    if relax
        optimize!(m)
        relaxation = JuMP.objective_value(m)
    end

    for i in 1:F
        for j in 1:F
            set_binary(xF[i, j])
        end
        for j in 1:R
            set_binary(xR[i, j])
        end
        set_binary(y[i])
        set_binary(xr[i])
    end

    set_attribute(m, "CPXPARAM_TimeLimit", 15 * 60)
    start_time = time()
    optimize!(m)
    end_time = time()

    if !relax
        relaxation = JuMP.objective_bound(m)
    end

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        return relaxation, JuMP.objective_value(m), JuMP.objective_bound(m), JuMP.node_count(m), end_time - start_time
    else
        return relaxation, false, JuMP.objective_bound(m), JuMP.node_count(m), end_time - start_time
    end
end

function runInstance(instanceName::String, p::Int)
    F, R, distanceMatrix = readInstance_tsp("MLA_Projet/tsp_data/" * instanceName * ".tsp")
    threshold = Int(floor(mean(distanceMatrix)))
    println("SCFF")
    relaxation, UB, LB, noeuds, t = SCFF(F, R, distanceMatrix, threshold, p, true, true, false)
    # relaxationDis, UBDis, LBDis, noeudsDis, tDis = SCFF(F, R, distanceMatrix, threshold, p, false, false, false)
    if UB == false
        gapLBUB = "_"
    else
        gapLBUB = (UB - LB)/UB
    end
    
    f = open("MLA_Projet/results/benchmark_grp9.csv", "a")
    write(f, instanceName)
    write(f, ";")
    write(f, string(floor(p)))
    write(f, ";")
    write(f, "SCFF")
    write(f, ";")
    if UB == false
        write(f, "_")
    else
        write(f, string(round(UB, digits=3)))
    end
    write(f, ";")
    write(f, string(round(LB, digits=3)))
    write(f, ";")
    if UB == false
        write(f, "_")
    else
        write(f, string(round(gapLBUB, digits=3)))
    end
    write(f, ";")
    write(f, string(round(t, digits=3)))
    write(f, ";")
    write(f, string(floor(noeuds)))
    write(f, ";")
    write(f, "_")
    write(f, ";")
    write(f, "_")
    write(f, "\n")
    close(f)

    println("MTZ")
    relaxation, UB, LB, noeuds, t = MTZ(F, R, distanceMatrix, threshold, p, true, true, false)
    if UB == false
        gapLBUB = "_"
    else
        gapLBUB = (UB - LB)/UB
    end
    
    f = open("MLA_Projet/results/benchmark_grp9.csv", "a")
    write(f, instanceName)
    write(f, ";")
    write(f, string(floor(p)))
    write(f, ";")
    write(f, "MTZ")
    write(f, ";")
    if UB == false
        write(f, "_")
    else
        write(f, string(round(UB, digits=3)))
    end
    write(f, ";")
    write(f, string(round(LB, digits=3)))
    write(f, ";")
    if UB == false
        write(f, "_")
    else
        write(f, string(round(gapLBUB, digits=3)))
    end
    write(f, ";")
    write(f, string(round(t, digits=3)))
    write(f, ";")
    write(f, string(floor(noeuds)))
    write(f, ";")
    write(f, "_")
    write(f, ";")
    write(f, "_")
    write(f, "\n")
    close(f)

end


F, R, distanceMatrix = readInstance_tsp("MLA_Projet/tsp_data/qa194.tsp")

threshold = Int(floor(mean(distanceMatrix)))
println("TEST")
println(threshold)
p = 10


#SCFF(F, R, distanceMatrix, threshold, p, true, true, false)

#MCFF(F, R, distanceMatrix, threshold, p, true, true, false)

MTZ(F, R, distanceMatrix, threshold, p, true, true, false)

# benchmarks = ["qa194", "lu980", "zi929"]
# for instance in benchmarks
#     for p in [10, 15, 20]
#         runInstance(instance, p)
#     end
# end
