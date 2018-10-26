# This script takes arguments from the command line which is assumed to be
# formatted in the way that
# julia driver.jl g ν H L T γ M Δt
# We also assume that the script is run from a work directory where we can put output
# files.


# Setup variables
#---------------------------------------------------------------------------------------------------

#const DT_MAX = 100000
const THERM_FRAC = 1/10                 # The fraction the thermalization time is multiplied by to
                                        # get time between measurements Δt
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
Δt = parse(Int64, ARGS[8])      # Re-interpreted as the maximum value Δt can have.

# Calculate remaining parameters

# Calculate periodic boundary conditioned f s.t. fL ∈ N
f = ceil(abs(H/(2π)*L))/L*sign(H)
# Updating H to value given by f
H = 2π*f
# Calculate inverse temperature
β = 1/T

# TODO: We should have some more checks and balances for the values entered in the variables

work_dir = "$(pwd())"
@everywhere println(pwd())
# Broadcast work_dir to all processes
@eval @everywhere work_dir=$work_dir
@everywhere source_dir = "$(work_dir)/../../Source"
println("driver.jl was called from $(pwd())")
println("Loading MC code")
@everywhere include("$(source_dir)/ChiralMC.jl")
@everywhere using ChiralMC

#@everywhere include("./functions_msc.jl")
@everywhere include("$(source_dir)/functions_observables.jl")
@everywhere include("$(source_dir)/functions_parallel.jl")
include("$(source_dir)/functions_plots_and_files.jl")
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


# Run simulation
#---------------------------------------------------------------------------------------------------

println("Running simulation..\n\nThermalizing states from high energy
--------------------------------------------------------------------------------------")

@time (t₀, ψ_ref, sim_ref, ψ_w, sim_w) = initializeParallelStatesS(syst, sim)
println("Initialization steps finished. Thermalization plots written to file.")

# Changed from original behavior: use t₀ to determine Δt
Δt = min(Δt, ceil(Int64, t₀*THERM_FRAC))
av_u⁺, av_u⁻ = meanAmplitudes(ψ_ref)
max_u⁺, min_u⁺, max_u⁻, min_u⁻ = maxMinAmplitudes(ψ_ref)
println("Amplitudes of reference state:\n⟨u⁺⟩ =\t$(av_u⁺)\t\t\t⟨u⁻⟩ =\t$(av_u⁻)\nmax(u⁺) =\t$(max_u⁺)\t\tmax(u⁻) =\t$(max_u⁻)
min(u⁺) =\t$(min_u⁺)\t\tmin(u⁻) =\t$(min_u⁻)")

# Setup files
#---------------------------------------------------------------------------------------------------

println("Writing system:\ng\t=\t$(g)\nν\t=\t$(ν)\nH\t=\t$(H)\nL\t=\t$(L)\nT\t=\t$(T)")
println("γ\t=\t$(γ)\nM\t=\t$(M)\nt₀\t=\t$(t₀)\nΔt\t=\t$(Δt)\nf\t=\t$(f)")
writeSimulationConstants(syst, sim, M, t₀, Δt)

# Saving the reference state and the first worker state
save(ψ_ref, "state_ref")
save(ψ_w[1], "state_w")
# Including ψ_ref in the list of states as the last state
ψ_list = vcat(ψ_w, [ψ_ref])


println("\nMeasuring structure function and vortex lattice
--------------------------------------------------------------------------------------")
@time (av_V⁺, err_V⁺, V⁺, av_V⁻, err_V⁻, V⁻, av_S⁺, err_S⁺, S⁺, av_S⁻, err_S⁻, S⁻) = parallelSFVLA!(k_matrix, ψ_list, sim_ref, M, Δt)
plotStructureFunctionVortexLatticeS(av_V⁺, av_V⁻, V⁺[rand(1:M)], V⁻[rand(1:M)], av_S⁺, av_S⁻, k_matrix)

println("\nSimulation finished!  o(〃＾▽＾〃)o\n\nResults are found in \n$(pwd())")
