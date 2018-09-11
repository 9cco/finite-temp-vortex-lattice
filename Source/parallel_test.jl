# This script takes arguments from the command line which is assumed to be
# formatted in the way that
# julia driver.jl g ν H L T γ M Δt
# We also assume that the script is run from the Source directory and that there is
# a neighboring Data directory.

tic()

# Setup variables
#---------------------------------------------------------------------------------------------------

# Go through argument list and check that we only have numbers
const variables = ["g", "ν", "H", "L", "T", "γ", "M", "Δt"]
@everywhere const two_pi = 2π

function helpMessage(variables::Array{String, 1})
	print("\ndriver.jl takes $(length(variables)) number of variables
use: julia driver.jl ")
	for i = 1:length(variables)
		print("$(variables[i]) ")
	end
	print("\n")
	println("Where each variable can be written e.g. in the format 1.32e-3\n")
end

length(ARGS) == length(variables) || throw(error("ERROR: Wrong number of arguments"))

# Scrubbing input
for i = 1:length(ARGS)
	if !ismatch(r"^[+-]?[0-9]*?[\.]?[0-9]*?[e]?[+-]?[0-9]*?$", ARGS[i])
		helpMessage(variables)
		println("(๑•̀ㅁ•́๑)✧
Input $(ARGS[i]) is an invalid number for $(variables[i])\n")
		throw(DomainError())
	end
end

# Parsing input as floats.
g = parse(Float64, ARGS[1])
ν = parse(Float64, ARGS[2])
H = parse(Float64, ARGS[3])
L = parse(Float64, ARGS[4])
T = parse(Float64, ARGS[5])
γ = parse(Float64, ARGS[6])
M = parse(Int64, ARGS[7])
Δt = parse(Int64, ARGS[8])

# Calculate remaining parameters

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = ceil(abs(H/(2π)*L))/L*sign(H)
# Calculate inverse temperature
β = 1/T

# TODO: We should have some more checks and balances for the values entered in the variables

println("Loading MC code")

@everywhere include("./ChiralMC.jl")
@everywhere using ChiralMC

#@everywhere include("./functions_msc.jl")
#@everywhere include("./functions_observables.jl")
@everywhere include("./functions_parallel.jl")
include("./functions_plots_and_files.jl")
# Hack to disable gr error messages.
ENV["GKSwstype"] = "100"

# Create system
syst = SystConstants(L, γ, 1/g^2, ν, f, β)
sim = Controls(π/3, 0.4, 3.0)

# Construct k-matrix where the horizontal axis contains kx ∈ [-π, π), while
# the vertical axis contain ky ∈ [-π, π) at the second component
k_matrix = [[2π/L*(x-1-L/2), 2π/L*(L/2-y)] for y=1:L, x=1:L];

println("\n\n
##########################################################################################

       Welcome to ChiralMC simulation software - your friend in superconductors

##########################################################################################
\nAll parameters seems to be in order, we are ready to begin\n
                       ｡･:*:･ﾟ★,｡･:*:･ﾟ☆ very d(*^▽^*)b good ｡･:*:･ﾟ★,｡･:*:･ﾟ☆\n\n")

# Setup files
#---------------------------------------------------------------------------------------------------

println("Setting up folder and system files")
# Make and create system directory in neighboring Data directory and move to this
mkcdSystemDirectory(syst, M, Δt)
writeSimulationConstants(syst, sim, M, Δt)

# Run test
#---------------------------------------------------------------------------------------------------

ψ = State(2, syst)
println("Starting test of $(M) parallel Monte-Carlo steps on $(nprocs()) processes")
ψ_list = parallelMultiplyState(ψ, sim, M)

println("\nSimulation finished!  o(〃＾▽＾〃)o\n\nResults are found in $(pwd())")
toc()
