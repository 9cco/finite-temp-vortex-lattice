# This script is made to make sure julia has updated all packages and
# installed the neccessary packages to run our simulation
#
using Pkg

# Print installed packages
println(Pkg.installed())
Pkg.resolve()

# Making sure necessary packages are included:
Pkg.add("Distributed")
Pkg.add("Test")
Pkg.add("BenchmarkTools")
Pkg.add("LinearAlgebra")
Pkg.add("LaTeXStrings")
Pkg.add("Dates")
Pkg.add("Primes")
Pkg.add("SharedArrays")
Pkg.add("DelimitedFiles")
Pkg.add("PyCall")
Pkg.add("PyPlot")
Pkg.add("NLsolve")
Pkg.add("JLD")
Pkg.add("StatsBase")
Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("MCMCDiagnostics")

# Updating packages
Pkg.update()
Pkg.build()

# Checking package status
println(Pkg.status())

using Plots
pyplot()
