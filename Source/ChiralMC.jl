module ChiralMC

using Distributions
using StatsBase
using Base.Test
using MCMCDiagnostics

export SystConstants, LatticeSite, State, Controls, NearestNeighbors, NextNeighbors, NNNeighbors, two_pi

include("types.jl")

####################################################################################################
#                                 Main section of functions
#
####################################################################################################

export maxRelErr, scientificRounding, avgErr, maxRelErrString, splitParallell, matrixIndices, findMaximaIndices

include("functions_msc.jl")

# This export is for testing
export latticeNeighbors, latticeNextNeighbors, latticeNNNeighbors
export checkNeighbors

include("functions_neighbors.jl")

export save, readState, checkState, set!, meanAmplitudes, maxMinAmplitudes, addToList, loadStates
include("functions_types.jl")

export E, ΔE

include("functions_energy.jl")

export mcSweep!, mcProposalFraction!, mcProposalFraction, mcSweepEn!, mcSweepFrac!, adjustSimConstants!, findEquilibrium
export nMCSEnergyDynamic, adjustSimConstantsPar, nMCSEnergy

include("functions_mc.jl")

export parallelThermalization!, findEquilibrium

include("functions_thermalization.jl")

end # Module ChiralMC
