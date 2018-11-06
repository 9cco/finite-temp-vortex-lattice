module ChiralMC

using Distributions
using StatsBase
using Base.Test

export SystConstants, LatticeSite, State, Controls, Neighbors, NearestNeighbors, NextNeighbors, NNNeighbors, two_pi

include("types.jl")

####################################################################################################
#                                 Main section of functions
#
####################################################################################################

include("functions_msc.jl")

# This export is for testing
export latticeNeighbors, latticeNextNeighbors, latticeNNNeighbors
export checkNeighbors

include("functions_neighbors.jl")

export save, checkState, set!, meanAmplitudes, maxMinAmplitudes
include("functions_types.jl")

export E, ΔE

include("functions_symmetric_energy.jl")

export mcSweep!, mcProposalFraction!, mcProposalFraction, mcSweepEn!, mcSweepFrac!, adjustSimConstants!, findEquilibrium
export nMCSEnergyDynamic, adjustSimConstantsPar

include("functions_mc.jl")

export parallelThermalization!, findEquilibrium

include("functions_thermalization.jl")

end # Module ChiralMC
