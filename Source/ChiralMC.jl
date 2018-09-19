module ChiralMC

using Distributions
using StatsBase
using Base.Test

export SystConstants, LatticeSite, State, Controls, Neighbors, NextNeighbors, NNNeighbors

include("types.jl")

####################################################################################################
#                                 Main section of functions
#
####################################################################################################

include("functions_msc.jl")

export save
include("functions_types.jl")

export E

include("functions_energy.jl")

export mcSweep!, mcProposalFraction!, mcSweepEn!, mcSweepFrac!, adjustSimConstants!, findEquilibrium

include("functions_mc.jl")

end # Module ChiralMC
