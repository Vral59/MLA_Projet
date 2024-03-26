using Pkg
Pkg.activate(".") # C'est mieux de pas polluer l'environnement de base
Pkg.add("JuMP")
Pkg.add("CPLEX")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("StatsBase")
Pkg.add("Statistics")
