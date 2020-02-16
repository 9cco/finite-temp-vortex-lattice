# Current intention:
# Using not as many thermalization steps like we figured out that we could use previously we want to
# explore a wide range of ν values at different weights of the MGT term.

# Change list: erland <-> fram
# out_path, staging_path, src_path, split, readline comment, target temp

# Change list: κs -> νs
# setting vector and parse argument, ARGS println, init_files, :%s/N_κ/N_ν/g, N_ν length, :%s/_by_κ/_by_ν/g, out_folder
# estimate AR, println MEASUREMENTS, save meta, save final state, :%s/κs/κ/gc, search manually for ν -> νs 

println("
  _________                                                   .___             __  .__      .__  __          
 /   _____/__ ________   ___________   ____  ____   ____    __| _/_ __   _____/  |_|__|__  _|__|/  |_ ___.__.
 \\_____  \\|  |  \\____ \\_/ __ \\_  __ \\_/ ___\\/  _ \\ /    \\  / __ |  |  \\_/ ___\\   __\\  \\  \\/ /  \\   __<   |  |
 /        \\  |  /  |_> >  ___/|  | \\/\\  \\__(  <_> )   |  \\/ /_/ |  |  /\\  \\___|  | |  |\\   /|  ||  |  \\___  |
/_______  /____/|   __/ \\___  >__|    \\___  >____/|___|  /\\____ |____/  \\___  >__| |__| \\_/ |__||__|  / ____|
        \\/      |__|        \\/            \\/           \\/      \\/           \\/                        \\/     
        ")
println("Script created by Fredrik Nicolai Krohg @ fredrik.n.krohg@ntnu.no")
println("This line updated on 18/12")

using Distributed
@everywhere using Distributions
using Test
#using BenchmarkTools
using Dates

@everywhere struct Hack end
function fixRC()
    for p in workers()
        @fetchfrom p Hack()
    end
end
fixRC()

@everywhere src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
#@everywhere src_path = "/cluster/home/fredrkro/mc/Source/Grid/"
out_path = "/home/nicolai/mc/Scripts/"
#out_path = "/cluster/home/fredrkro/mc/Data/FieldRuns/"
staging_path = ""
#staging_path = "/cluster/home/fredrkro/StagingField/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")
include("script_functions.jl")

using JLD

# Enter data directory for structure function
fixRC()
# We run a simulation with parameter ν supplied by script
if length(ARGS) != 1
    println("ERROR: Need κ supplied as script argument.")
    exit(-1)
end
κ = parse(Float64, ARGS[1])                         # MGT weights
νs =  [0.3] #[-0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.4, 0.5, 0.7]
g = 0.3                               # Gauge couplings
κ₅ = 1.0
N_ν = length(νs)

# Other parameters
M_est = 2^6  # Number of MC-steps to do at T_start before cooldown.
M_col = 2^16 # Number of MC-steps to use for cooling down the systems
N_steps = 2^10   # Number of temperatures to go through from high temp before reaching the final temperature. (will divide M_col) must be >= 2
M_th = 2^16  # Number of MC-steps to do for thermalization.
M = 2^13     # Number of measurements
M_amp = 2^3   # Number of these measurements that will be amplitude measurements, i.e. we need M >= M_amp
Δt = 2^6      # Number of MC-steps between each measurement
# L is assumed to be even.
L = 42     # System length
L₁ = L
L₂ = L
L₃ = L
split = (1,1,4) #(1,1,2)#
N = L₁*L₂*L₃
# Specify extended landau gauge using n,m s.t. constant Aᵢ(Lⱼ+1) = Aᵢ(0) mod 2π
# n | L₁ and m | L₂
n = 1; m = 0
f = n/L₁ - m/L₂
Ts = [1.756, 1.752, 1.748, 1.744, 1.741, 1.736, 1.732, 1.728, 1.724, 1.721, 1.716, 1.712, 1.708, 1.704, 1.701]#[1.805, 1.801, 1.795, 1.79, 1.785, 1.781, 1.775, 1.77, 1.765, 1.761]#[1.845, 1.84, 1.835, 1.83, 1.825, 1.82, 1.815, 1.81]#[1.45, 1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.95, 0.9] #[4.0, 2.5, 2.0, 1.9, 1.8, 1.7, 1.65, 1.6, 1.55, 1.5]#[2.0, 1.9, 1.8, 1.7]
T_start = 2*4.0+0.1   # Start-temperature of states when thermalizing from ab-inito state


# Print parameters
println("\nStarted script at: $(Dates.format(now(), "HH:MM")) on $(Dates.format(now(), "d/m"))")
println("Target temps: $(Ts)")
println("f = $(f)")
println("n = $(n)")
println("m = $(m)")
println("L = $(L)")
println("νs = $(νs)")
println("κ₅ = $(κ₅)")
println("κ = $(κ)")
println("g = $(g)")
println("Each state split in $(split)")
if f < 0
    println("Number of anti-vortices: $(-f*L₁*L₂)")
else
    println("Number of vortices: $(f*L₁*L₂)")
end
flush(stdout)
cd(out_path)

# Checking available workers
workers_pr_state = Int(prod(split))
needed_workers = workers_pr_state*N_ν - 1 + 1# + N_ν
if nprocs()-1 < needed_workers
    #println("ERROR: Only $(nprocs()-1) workers have been allocated. $(needed_workers) are needed to split all
#of the $((needed_workers-N_ν+1)/workers_pr_state) states into $(workers_pr_state) subcuboids and then
#run these in parallel using $(N_ν) additional workers.")
    println("ERROR: Only $(nprocs()-1) workers have been allocated. $(needed_workers) are needed to split all
of the $((needed_workers-1+1)/workers_pr_state) states into $(workers_pr_state) subcuboids.")
    exit()
end



# Initiating sates
################################################################################################
init_files = [staging_path*"init_state_g=$(round(g; digits=3))_nu=$(round(ν; digits=3))_kap=$(round(κ, digits=3)).jld" for ν in νs]
cubs, systs, init_temps = initiateStates(init_files, L₁, L₂, L₃, g, νs, κ₅, κ, n, m, split)



# Time-estimation
################################################################################################

N_T = length(Ts)
t_MCS = estimateRuntime!(cubs, N, M_est, Δt)

# Estimating ETC
println("Measurements and thermalization will do $(N_T*(Δt*M + M_th + M_col)) MCS for each of the $(N_ν)
states, which has an ETC of $(round((Δt*M + M_th + M_col)*N_T*N_ν*t_MCS/3600; digits=2)) h")
print("Continue with thermalization and measurements? (y/n): ")
flush(stdout)

#user_input = readline(stdin)
#if user_input == "n"
#    exit()
#end

for T in Ts         #>>>>>>>>>>>>>>>>> Start of T loop <<<<<<<<<<<<<<<<<<<<#
# Making global variables be assignable in the new local scope
global M_col

out_folders = ["C4_model_L=$(L)_T=$(round(T; digits=3))_g=$(round(g; digits=3))_nu=$(round(ν; digits=3))_kap=$(κ)_n=$(n)_m=$(m)_mult_kap" for ν in νs]
for folder in out_folders
    mkcd(folder)
    cd("../")
    println("Making folder $(folder) in path $(out_path)")
end


println("\nCOOLDOWN & THERMALIZATION T -> $(T)\n#############################################################")

# Cooldown
################################################################################################

M_col, t_col = cooldownStates!(cubs, N, M_col, N_steps, M_est, init_temps, T, out_folders)
# Setting initial temps to be the current T so that next time we go through the T for loop we will
# cool down from the previous temperature.
fill!(init_temps, T)

t_MCS0 = t_col/(M_col*N_ν)



# Thermalization
################################################################################################

t_MCS1 = thermalizeStates!(cubs, M_th, M_est, M_col, out_folders)

t_MCS_av = (t_MCS0 + t_MCS1)/2

AR_est = [estimateAR!(cub)[1] for cub in cubs]
println("After thermalization the estimated AR are:")
for (k, est) = enumerate(AR_est); println("AR(ν=$(round(νs[k]; digits=1))): ≈ $(round(est*100; digits=2))%"); end



# Doing measurements
################################################################################################
println("\nMEASUREMENTS\n#############################################################")
println("We will now do $(Δt*M) MCS pr ν value  which gives
$(M) measurements and has an ETC of $(round(Δt*M*N_ν*t_MCS_av/3600; digits=2)) h")
println("Started measurements at: $(Dates.format(now(), "HH:MM")) on $(Dates.format(now(), "d/m"))")
flush(stdout)

# Setup storage for dual stifness and energy
E_by_ν = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
E_prev = [energy(cub) for cub in cubs]
# Storage for vortices and dual vortices
S⁺_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
S⁻_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
V⁺_by_ν = [Array{Array{Float64, 2}, 1}(undef, M) for k = 1:N_ν]
V⁻_by_ν = [Array{Array{Float64, 2}, 1}(undef, M) for k = 1:N_ν]
Vx_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
Vy_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
Sx_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
Sy_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
# Storage for amplitude matrices
u⁺_lattices_by_ν = [Array{Array{Float64, 3},1}(undef, M_amp) for k = 1:N_ν]
u⁻_lattices_by_ν = [Array{Array{Float64, 3},1}(undef, M_amp) for k = 1:N_ν]
u⁺_xy_lattices_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
u⁻_xy_lattices_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
u⁺_avg_by_ν = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
u⁻_avg_by_ν = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
δu²_by_ν = [Array{Float64, 1}(undef, M) for k = 1:N_ν] # 1/N*sum_r( (u⁺ᵣ)^2 - (u⁻ᵣ)^2 )

Bz_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]
Δθ_avg_by_ν = [zeros(Float64, L₁, L₂) for k = 1:N_ν]

# Storage for helicity modulus derivatives
#dH_01_x = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_01_y = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_01_z = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_10_x = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_10_y = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_10_z = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_11_x = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_11_y = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#dH_11_z = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_01_x = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_01_y = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_01_z = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_10_x = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_10_y = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_10_z = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_11_x = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_11_y = [Array{Float64, 1}(undef, M) for k = 1:N_ν]
#d²H_11_z = [Array{Float64, 1}(undef, M) for k = 1:N_ν]


# We preform M measurements
t_meas = @elapsed @time for m = 1:M
    dE_lists = nMCSEnUp!(cubs, Δt)[1]
    δE = [sum(dE_list) for dE_list in dE_lists]

    # Sort results in terms of temperature.
    for i = 1:N_ν
        E_by_ν[i][m] = E_prev[i] + δE[i]
        E_prev[i] += δE[i]
    end

    if m <= M_amp
        (proj_V⁺, proj_V⁻, S⁺, S⁻, proj_Vx, proj_Vy, Sx, Sy, u⁺_lattices, u⁻_lattices, u⁺_avg_z, u⁻_avg_z, u⁺_avg, u⁻_avg,
         δu², Bz_proj, Δθ_proj) = measureStates(cubs; full_amplitude_lattice = true)
#         dH_01s, dH_10s, dH_11s, d²H_01s, d²H_10s, d²H_11s) = measureStates(cubs; full_amplitude_lattice = true)
    else
        (proj_V⁺, proj_V⁻, S⁺, S⁻, proj_Vx, proj_Vy, Sx, Sy, u⁺_avg_z, u⁻_avg_z, u⁺_avg, u⁻_avg,
         δu², Bz_proj, Δθ_proj) = measureStates(cubs, full_amplitude_lattice = false)
#         dH_01s, dH_10s, dH_11s, d²H_01s, d²H_10s, d²H_11s) = measureStates(cubs, full_amplitude_lattice = false)
    end

    # Organize measurements into storage arrays
    for k = 1:N_ν
        V⁺_by_ν[k][m] = proj_V⁺[k]; V⁻_by_ν[k][m] = proj_V⁻[k]
        S⁺_avg_by_ν[k] .+= S⁺[k]; S⁻_avg_by_ν[k] .+= S⁻[k]
        Vx_avg_by_ν[k] .+= proj_Vx[k]; Vy_avg_by_ν[k] .+= proj_Vy[k]
        Sx_avg_by_ν[k] .+= Sx[k]; Sy_avg_by_ν[k] .+= Sy[k]

        # Save z-averaged lattices
        u⁺_xy_lattices_by_ν[k] .+= u⁺_avg_z[k]
        u⁻_xy_lattices_by_ν[k] .+= u⁻_avg_z[k]
        u⁺_avg_by_ν[k][m] = u⁺_avg[k]
        u⁻_avg_by_ν[k][m] = u⁻_avg[k]
        if m <= M_amp
            u⁺_lattices_by_ν[k][m] = u⁺_lattices[k]
            u⁻_lattices_by_ν[k][m] = u⁻_lattices[k]
        end
        δu²_by_ν[k][m] = δu²[k]

        Bz_avg_by_ν[k] .+= Bz_proj[k]
        Δθ_avg_by_ν[k] .+= Δθ_proj[k]

#        dH_01_x[k][m], dH_01_y[k][m], dH_01_z[k][m] = dH_01s[k]
#        dH_10_x[k][m], dH_10_y[k][m], dH_10_z[k][m] = dH_10s[k] 
#        dH_11_x[k][m], dH_11_y[k][m], dH_11_z[k][m] = dH_11s[k]
#
#        d²H_01_x[k][m], d²H_01_y[k][m], d²H_01_z[k][m] = d²H_01s[k] 
#        d²H_10_x[k][m], d²H_10_y[k][m], d²H_10_z[k][m] = d²H_10s[k]
#        d²H_11_x[k][m], d²H_11_y[k][m], d²H_11_z[k][m] = d²H_11s[k]
    end
end
S⁺_avg_by_ν = S⁺_avg_by_ν./M; S⁻_avg_by_ν = S⁻_avg_by_ν./M
Vx_avg_by_ν = Vx_avg_by_ν./M; Vy_avg_by_ν = Vy_avg_by_ν./M
Sx_avg_by_ν = Sx_avg_by_ν./M; Sy_avg_by_ν = Sy_avg_by_ν./M
Bz_avg_bu_ν = Bz_avg_by_ν./M; Δθ_avg_by_ν = Δθ_avg_by_ν./M
final_lattices = [getLattice(cub) for cub in cubs]
vortices_by_ν = [vortexSnapshot(lattice, getSyst(cubs[1])) for lattice in final_lattices]

println("Measurements used $(round(t_meas/3600; digits=1)) h, which means it took on average $(round(t_meas/(M*N_ν); digits=1)) s pr. measurement.")
println("and $(round(t_meas/(M*Δt*N*N_ν)*1e6; digits=2)) μs pr. MCS pr. lattice site.")
println("-------------------------------------------------------------\n\n")



# Saving results to file
################################################################################################
for k = 1:N_ν
    JLD.save(out_folders[k]*"/energies.jld", "Es", E_by_ν[k])

    JLD.save(out_folders[k]*"/meta.jld", "L1", L₁, "L2", L₂, "L3", L₃, "M", M, "M_amp", M_amp, "dt", Δt, "T", T, "n", n, "m", m, "kap5", κ₅, "kap", κ, "nu", νs[k], "g", g)
    JLD.save(out_folders[k]*"/vorticity.jld", "vortexes", vortices_by_ν[k], "sp_avg", S⁺_avg_by_ν[k], "sm_avg", S⁻_avg_by_ν[k], "vp", V⁺_by_ν[k], "vm", V⁻_by_ν[k])
    JLD.save(out_folders[k]*"/XY_vorticity.jld", "vx_avg", Vx_avg_by_ν[k], "vy_avg", Vy_avg_by_ν[k], "sx_avg", Sx_avg_by_ν[k], "sy_avg", Sy_avg_by_ν[k])
    JLD.save(out_folders[k]*"/amplitudes.jld", "up_lattices", u⁺_lattices_by_ν[k], "um_lattices", u⁻_lattices_by_ν[k], "up_xy", u⁺_xy_lattices_by_ν[k], "um_xy", u⁻_xy_lattices_by_ν[k], "up_avg", u⁺_avg_by_ν[k], "um_avg", u⁻_avg_by_ν[k], "du2", δu²_by_ν[k])
    JLD.save(out_folders[k]*"/vortex_consequences.jld", "Bz_avg", Bz_avg_by_ν[k], "phaseDiff_avg", Δθ_avg_by_ν[k])

    # Saving Helicity moduli
#    JLD.save(out_folders[k]*"/hel_mod.jld", "dH_01_x", dH_01_x[k], "dH_01_y", dH_01_y[k], "dH_01_z", dH_01_z[k],
#             "dH_10_x", dH_10_x[k], "dH_10_y", dH_10_y[k], "dH_10_z", dH_10_z[k],
#             "dH_11_x", dH_11_x[k], "dH_11_y", dH_11_y[k], "dH_11_z", dH_11_z[k],
#             "d2H_01_x", d²H_01_x[k], "d2H_01_y", d²H_01_y[k], "d2H_01_z", d²H_01_z[k],
#             "d2H_10_x", d²H_10_x[k], "d2H_10_y", d²H_10_y[k], "d2H_10_z", d²H_10_z[k],
#             "d2H_11_x", d²H_11_x[k], "d2H_11_y", d²H_11_y[k], "d2H_11_z", d²H_11_z[k])

    # Saving final state
    JLD.save(out_folders[k]*"/final_state_g=$(round(g; digits=3))_nu=$(round(νs[k]; digits=3))_kap=$(round(κ, digits=3)).jld", "lattice", final_lattices[k], "syst", getSyst(cubs[k]), "T", getTemp(cubs[k]), "controls", getControls(cubs[k]))
end

end          #>>>>>>>>>>>>>>>>> End of T loop <<<<<<<<<<<<<<<<<<<<#

println("ヾ(＾-＾)ノ")
