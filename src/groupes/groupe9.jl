# Salma Hammani et Francesco Gallo
# Problème du p-centre connexe
using JuMP
using CPLEX
using Statistics

include("../readData.jl")

function SCFF(F::Int, R::Int, distanceMatrix::Matrix{Int}, threshold::Int, p::Int, verbose=true)
    m = Model(CPLEX.Optimizer)
    if !verbose
        set_silent(m)
    end
    @variable(m, xF[1:F, 1:F], Bin) # x^F_ij avec i dans F et j and F
    @variable(m, xr[1:F], Bin) # x^F_ri avec i dans F et r la racine
    @variable(m, y[1:F], Bin) # y_i avec i dans F
    @variable(m, yr, Bin) # y_r où r est la racine
    @variable(m, xR[1:F, 1:R], Bin) # x^R_ij avec i dans F et j dans R
    @variable(m, r >= 0) # le rayon maximal
    @variable(m, f[1:F, 1:F] >= 0) # f_ij avec i et j dans F
    @variable(m, fr[1:F] >= 0) # f_ri avec i dans F

    @constraint(m, [j in 1:R], sum(xR[i, j] for i in 1:F) == 1)
    @constraint(m, [j in 1:R, i in 1:F], xR[i, j] <= y[i])
    @constraint(m, sum(y[i] for i in 1:F) == p)
    @constraint(m, [j in 1:R], r >= sum(distanceMatrix[i, j] * xR[i, j] for i in 1:F))

    @constraint(m, yr == 1)
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

    @objective(m, Min, r)

    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        println("Valeur de l'objectif final : ", JuMP.objective_value(m))
    else
        println("Pas de solution possible")
    end
end

function MCFF(F::Int, R::Int, distanceMatrix::Matrix{Int}, threshold::Int, p::Int, verbose=true)
    m = Model(CPLEX.Optimizer)
    if !verbose
        set_silent(m)
    end
    @variable(m, xF[1:F, 1:F], Bin) # x^F_ij avec i dans F et j and F
    @variable(m, xr[1:F], Bin) # x^F_ri avec i dans F et r la racine
    @variable(m, y[1:F], Bin) # y_i avec i dans F
    @variable(m, yr, Bin) # y_r où r est la racine
    @variable(m, xR[1:F, 1:R], Bin) # x^R_ij avec i dans F et j dans R
    @variable(m, r >= 0) # le rayon maximal
    @variable(m, fk[1:F, 1:F,1:F] >= 0) # f_k_ij avec i et j dans F
    @variable(m, frk[1:F,1:F] >= 0) # fk_ri avec i dans F


    @constraint(m, [j in 1:R], sum(xR[i, j] for i in 1:F) == 1)
    @constraint(m, [j in 1:R, i in 1:F], xR[i, j] <= y[i])
    @constraint(m, sum(y[i] for i in 1:F) == p)
    @constraint(m, [j in 1:R], r >= sum(distanceMatrix[i, j] * xR[i, j] for i in 1:F))

    @constraint(m, yr == 1)
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

    @objective(m, Min, r)

    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        println("Valeur de l'objectif final : ", JuMP.objective_value(m))
    else
        println("Pas de solution possible")
    end
end


function MTZ(F::Int, R::Int, distanceMatrix::Matrix{Int}, threshold::Int, p::Int, verbose=true)
    m = Model(CPLEX.Optimizer)
    if !verbose
        set_silent(m)
    end
    @variable(m, xF[1:F, 1:F], Bin) # x^F_ij avec i dans F et j and F
    @variable(m, xr[1:F], Bin) # x^F_ri avec i dans F et r la racine
    @variable(m, y[1:F], Bin) # y_i avec i dans F
    @variable(m, yr, Bin) # y_r où r est la racine
    @variable(m, xR[1:F, 1:R], Bin) # x^R_ij avec i dans F et j dans R
    @variable(m, r >= 0) # le rayon maximal
    @variable(m, u[1:F] >= 0) 
    @variable(m, ur >= 0) 

    @constraint(m, [j in 1:R], sum(xR[i, j] for i in 1:F) == 1)
    @constraint(m, [j in 1:R, i in 1:F], xR[i, j] <= y[i])
    @constraint(m, sum(y[i] for i in 1:F) == p)
    @constraint(m, [j in 1:R], r >= sum(distanceMatrix[i, j] * xR[i, j] for i in 1:F))

    @constraint(m, yr == 1)
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
    @constraint(m, [i in 1:F,j in 1:F], (F+1) * xF[i,j] + u[i] <= u[j] + F )
    @constraint(m, [j in 1:F], (F+1) * xr[j] + ur <= u[j] + F )
    @constraint(m, sum(xr[i] for i in 1:F) == 1)

    @objective(m, Min, r)

    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        println("Valeur de l'objectif final : ", JuMP.objective_value(m))
    else
        println("Pas de solution possible")
    end
end


F, R, distanceMatrix = readInstance_tsp("MLA_Projet/tsp_data/wi29.tsp")

threshold = Int(floor(mean(distanceMatrix)))
println(threshold)
p = 5

#SCFF(F, R, distanceMatrix, threshold, p)

#MCFF(F, R, distanceMatrix, threshold, p)

MTZ(F, R, distanceMatrix, threshold, p)
