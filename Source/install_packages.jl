# This script is made to make sure julia has updated all packages and
# installed the neccessary packages to run our simulation

# Print installed packages
println(Pkg.installed())
Pkg.resolve()

# Making sure necessary packages are included:
Pkg.add("StatsBase")
Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("MCMCDiagnostics")

# Updating packages
Pkg.update()

# Checking package status
Pkg.status()
