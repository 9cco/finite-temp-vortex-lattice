# Current intention:
# Run symmetric extended Landau gauge to check if we still get striped patterns
# close to Tc for otherwise same parameters as previously.

# Change list: erland <-> fram
# out_path, staging_path, src_path, split, readline comment, target temp

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

#@everywhere src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
@everywhere src_path = "/cluster/home/fredrkro/mc/Source/Grid/"
#out_path = "/home/nicolai/mc/Scripts/"
out_path = "/cluster/home/fredrkro/mc/Data/"
#staging_path = ""
staging_path = "/cluster/home/fredrkro/Staging/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")

using JLD

# Enter data directory for structure function
fixRC()
# We run a simulation with parameter ν supplied by script
if length(ARGS) != 1
    println("ERROR: Need g supplied as script argument.")
    exit(-1)
end
ν = parse(Float64, ARGS[1])         # Gauge coupling
gs = [0.3]#[-0.5, 0.0, 0.3, 0.5, 0.7]         # Fermi surface anisotropy
κ₅ = 1.0

# TODO: We should probably benchmark vortexSnapshot and see that it is faster than mcSweeps
# Other parameters
M_est = 2^7  # Number of MC-steps to do at T_start before cooldown.
M_col = 2^18 # Number of MC-steps to use for cooling down the systems
N_steps = 2^10   # Number of temperatures to go through from high temp before reaching the final temperature. (will divide M_col) must be >= 2
M_th = 2^18  # Number of MC-steps to do for thermalization.
M = 2^12     # Number of measurements
M_amp = 30   # Number of these measurements that will be amplitude measurements, i.e. we need M >= M_amp
Δt = 2^6      # Number of MC-steps between each measurement
# L is assumed to be even.
L = 2*32     # System length
L₁ = L
L₂ = L
L₃ = L
split = (2,4,4)
N = L₁*L₂*L₃
# Specify extended landau gauge using n,m s.t. constant Aᵢ(Lⱼ+1) = Aᵢ(0) mod 2π
# n | L₁ and m | L₂
n = 1; m = -1
f = n/L₁ - m/L₂
T = 2.0
T_start = 2*4.0+0.1   # Start-temperature of states when thermalizing. 
N_g = length(gs)


# Print parameters
println("\nStarted script at: $(Dates.format(now(), "HH:MM"))")
println("Target temp: $(T)")
println("f = $(f)")
println("n = $(n)")
println("m = $(m)")
println("L = $(L)")
println("ν = $(ν)")
println("κ₅ = $(κ₅)")
println("gs = $(gs)")
if f < 0
    println("Number of anti-vortices: $(-f*L₁*L₂)")
else
    println("Number of vortices: $(f*L₁*L₂)")
end
flush(stdout)
cd(out_path)
out_folders = ["full_model_L=$(L)_T=$(round(T; digits=2))_g=$(round(g; digits=3))_nu=$(round(ν; digits=3))_n=$(n)_m=$(m)_mult_g" for g in gs]

# Checking available workers
workers_pr_state = Int(split[1]*split[2]*split[3])
needed_workers = workers_pr_state*N_g
if nprocs()-1 < needed_workers
    println("ERROR: Only $(nprocs()-1) workers have been allocated. $(needed_workers) are needed to split all
of the $(length(gs)) states into $(split[1]*split[2]*split[3]) subcuboids.")
    exit()
end



# Initiating states
################################################################################################
println("\nINITIATING STATES\n#############################################################")
#
# We check if the staging folder (could be current folder) includes an initial state file.
# We assume that initial states for all gs are included in this file.
init_file = staging_path*"init_state_nu=$(round(ν; digits=3))_L=$(L).jld"
if isfile(init_file)
    # Reading the initial file.
    init_file_di = JLD.load(init_file)
    init_lattices = init_file_di["lattices"]
    init_T = init_file_di["T"]
    init_n = init_file_di["n"]
    init_m = init_file_di["m"]
    init_κ₅ = init_file_di["kappa"]
    init_ν = init_file_di["nu"]
    init_gs = init_file_di["gs"]
    init_controls = init_file_di["controls"]
    init_lattice = init_lattices[1]

    N_g_init = length(init_gs)
    if (init_gs != gs || init_n != n || init_m != m || init_ν != ν || init_κ₅ != κ₅
        || L₁ != size(init_lattice, 1) || L₂ != size(init_lattice, 2) || L₃ != size(init_lattice, 3))
        println("ERROR: Parameters in initial-states file does not match target parameters in script.")
        exit(-1)
    end

    # Construct initial states
    cubs = Array{Cuboid, 1}(undef, N_g)
    systs = [SystConstants(L₁,L₂,L₃,1/init_g^2,init_ν,init_κ₅,init_n,init_m) for init_g in init_gs]
    init_temps = [init_T for syst in systs]
    for i = 1:N_g
        lattice = init_lattices[i]
        cubs[i] = Cuboid(lattice, systs[i], init_controls[i], split, 1/init_T; pid_start=(i-1)*workers_pr_state + N_g + 1)
    end
    
    println("$(N_g) state(s) created from initiation file with g-value(s):")
    println(init_gs)
    println("Controls of read state:")
    printControls(cubs)
    
    # Tuning initial controls
    println("Adjusting simulation constants in initial states.")
    flush(stdout)
    for k = 1:N_g; setUpdates!(cubs[k]; θ_max = 3.13, u_max = 1.0/√(2), A_max = 1e-1); end
    #@time for k = 1:N_g; retuneUpdates!(cubs[k]; θ_max = 0.2, u_max = 1.0/√(2), A_max = 9e-3); end
    println("Controls after adjustment:")
    printControls(cubs)
    AR_est = [estimateAR!(cub)[1] for cub in cubs]
    println("AR with adjusted controls:")
    for est in AR_est; println("AR: ≈ $(round(est*100; digits=2))%"); end

else # If no initial file is found, then we contruct initial states at temperature T_start and thermalize them
    # Make ab inito un-correlated phases state if no initial_state is given.
    systs = [SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,n,m) for g in gs]
    cubs = Array{Cuboid, 1}(undef, N_g)
    for k = 1:N_g
        cubs[k] = Cuboid(6, systs[k], split, 1/T_start; u⁺=1/√(2), u⁻=1/√(2), pid_start=(k-1)*workers_pr_state + N_g + 1)
    end
    init_temps = [T_start for syst in systs]

    println("No initialization file found. Thermalizing ab-inito states.")
    flush(stdout)
    # Tune and thermalize start-temperature
    t_init_th = @elapsed nMCS!(cubs, ceil(Int64, M_th/4))

    println("Adjusting simulation constants in initial states.")
    #@time for cub in cubs; tuneUpdates!(cub); end
    for k = 1:N_g; setUpdates!(cubs[k]; θ_max = 3.13, u_max = 1.0/√(2), A_max = 1e-1); end # If we set updates manually
    println("Controls after adjustment:")
    printControls(cubs)
    flush(stdout)

    t_init_th2 = @elapsed nMCS!(cubs, ceil(Int64, 3*M_th/4))
    println("Thermalization of ab-inito states used $(round((t_init_th+t_init_th2)/3600; digits=1)) h.")
end

for folder in out_folders
    mkcd(folder)
    cd("../")
    println("Making folder $(folder) in path $(out_path)")
end
println("\nCOOLDOWN & THERMALIZATION\n#############################################################")


# Time-estimation
################################################################################################

# Do M_est MC-sweeps on first
println("Estimating runtime")
t_meas = @elapsed for i = 1:M_est
    nMCSEnUp!(cubs, Δt)

    all_res = [dualStiffnesses(cub) for cub in cubs]
    En_res = [energy(cub) for cub in cubs]
    for cub in cubs; chiralAmplitudeSnapshot(cub); xyVortexSnapshot(cub); end
end
t_meas = t_meas/M_est
t_MCS = t_meas/(Δt*length(cubs))
println("We did $(M_est) measurement steps using average $(round(t_meas; digits=1)) s pr measurement. 
This yields $(round(t_MCS*1000; digits=1)) ms on average for each MCS")
println("and $(round(t_MCS/N*1e6; digits=1)) μs on average pr. lattice site.")

# Estimating ETC
println("Measurements and thermalization will do $(Δt*M + M_th + M_col) MCS for each of the $(N_g)
g-values, which has an ETC of $(round((Δt*M + M_th + M_col)*N_g*t_MCS/3600; digits=2)) h")
print("Continue with thermalization and measurements? (y/n): ")
flush(stdout)

#user_input = readline(stdin)
#if user_input == "n"
#    exit()
#end



# Cooldown
################################################################################################
println("\nStarted cooldown at: $(Dates.format(now(), "HH:MM"))")

# Do M_th MCS in N_steps steps
M_pr_step = ceil(Int64, M_col/N_steps)
M_col = M_pr_step*N_steps
# Each step is associated with a different set of temperatures given by each row in the matrix
temp_mt = genGeometricTemperatureSteps(init_temps, fill(T, N_g), N_steps)
# Energy storage
E_matrix = Array{Float64, 2}(undef, M_col+1, N_g);
for (i, cub) = enumerate(cubs); E_matrix[1,i] = energy(cub); end
E_therm = Array{Float64, 2}(undef, M_th, N_g);
# Vorticity storage
S⁺_col_by_g_aux = [Array{Array{Float64, 2}, 1}(undef, M_pr_step) for i = 1:N_g]
S⁻_col_by_g_aux = [Array{Array{Float64, 2}, 1}(undef, M_pr_step) for i = 1:N_g]
S⁺_col_by_g = [Array{Array{Float64, 2}, 1}(undef, N_steps) for i = 1:N_g]
S⁻_col_by_g = [Array{Array{Float64, 2}, 1}(undef, N_steps) for i = 1:N_g]
# Acceptance rate storage
accepts_matrix = Array{Int64, 2}(undef, M_col, N_g);

println("Cooling down from T=$(init_temps) using $(N_steps) steps with $(M_pr_step) MCS pr. step.")
flush(stdout)

t_col = @elapsed for step = 1:N_steps
    # First set the temperatures associated with this step to each system.
    for (k, cub) = enumerate(cubs); setTemp!(cub, temp_mt[step, k]); end

    # Do the MCS associated with this step
    for i = 1:M_pr_step
        δE_lists, accepts_lists = nMCSEnUp!(cubs, 1)
        for (k, cub) = enumerate(cubs)
            δE = δE_lists[k][1]
            accepts = accepts_lists[k][1]
            δE, accepts = mcSweepEnUp!(cub)
            E_matrix[(step-1)*M_pr_step+i+1, k] = E_matrix[(step-1)*M_pr_step+i, k] + δE
            accepts_matrix[(step-1)*M_pr_step+i, k] = accepts
            proj_V⁺, proj_V⁻ = xyVortexSnapshot(cub)
            S⁺_col_by_g_aux[k][i], S⁻_col_by_g_aux[k][i] = structureFunction(proj_V⁺, proj_V⁻)
        end
    end

    # In the end of a step we correct energy rounding-error
    for (j, cub) = enumerate(cubs); E_matrix[step*M_pr_step+1, j] = energy(cub); end

    # Then we also calculate the average structure function at this step.
    for k = 1:N_g;  S⁺_col_by_g[k][step] = mean(S⁺_col_by_g_aux[k]); S⁻_col_by_g[k][step] = mean(S⁻_col_by_g_aux[k]); end

end

println("Cooldown used $(round(t_col/3600; digits=1)) h, i.e. $(round(t_col/(M_col*N_g*N)*1e6; digits=2)) μs pr. lattice site")
E_matrix = E_matrix./N;

# Saving Cooldown to separate ν-folders
for k = 1:N_g
    JLD.save(out_folders[k]*"/cooldown.jld", "sp", S⁺_col_by_g[k], "sm", S⁻_col_by_g[k], "E_list", E_matrix[:, k], "accepts_list", accepts_matrix[:,k], "temp_list", temp_mt[:, k], "N_steps", N_steps, "M_pr_step", M_pr_step, "M_est", M_est)
end

t_MCS0 = t_col/(M_col*N_g)



# Adjusting simulation constants
################################################################################################

#println("Adjusting simulation constants.")
#@time for cub in cubs; tuneUpdates!(cub); end
#println("Controls after adjustment:")
#printControls(cubs)



# Thermalization
################################################################################################
println("Thermalizing at target temperatures for additional $(M_th) MCS")
println("\nStarted thermalization at: $(Dates.format(now(), "HH:MM"))")
flush(stdout)

for j = 1:N_g; E_therm[1, j] = energy(cubs[j]); end
t_th = @elapsed δE_lists = nMCSEnUp!(cubs, M_th-1)[1]
for i = 2:M_th
    for k = 1:N_g
        E_therm[i, k] = E_therm[i-1, k] + δE_lists[k][i-1]
    end
end
AR_est = [estimateAR!(cub)[1] for cub in cubs]

E_therm = E_therm./N

# Update time-estimation
t_MCS1 = t_th/(M_th*N_g)
t_MCS_av = (t_MCS0 + t_MCS1)/2

therm_lattices = [getLattice(cub) for cub in cubs]

# Saving thermalization
for k = 1:N_g
    JLD.save(out_folders[k]*"/thermalization.jld", "e_thm", E_therm[:,k], "syst", systs[k], "lattice", therm_lattices[k], "M_th", M_th, "M_est", M_est, "M_col", M_col)
end
println("Thermalization used $(round(t_th/3600; digits=1)) h, i.e. $(round(t_MCS1/(N_g*N)*1e6; digits=2)) μs pr. lattice site. Energies and states saved to file.")
println("After thermalization the estimated AR are:")
for (k, est) = enumerate(AR_est); println("AR(g=$(round(gs[k]; digits=1))): ≈ $(round(est*100; digits=2))%"); end



# Doing measurements
################################################################################################
println("\nMEASUREMENTS\n#############################################################")
println("We will now do $(Δt*M) MCS pr ν value  which gives
$(M) measurements and has an ETC of $(round(Δt*M*N_g*t_MCS_av/3600; digits=2)) h")
println("Started measurements at: $(Dates.format(now(), "HH:MM"))")
flush(stdout)

# Setup storage for dual stifness and energy
E_by_g = [Array{Float64, 1}(undef, M) for k = 1:N_g]
E_prev = [energy(cub) for cub in cubs]
# Storage for vortices and dual vortices
S⁺_by_g = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_g]
S⁻_by_g = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_g]
V⁺_avg_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
V⁻_avg_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
Vx_avg_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
Vy_avg_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
Sx_avg_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
Sy_avg_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
# Storage for amplitude matrices
u⁺_lattices_by_g = [Array{Array{Float64, 3},1}(undef, M_amp) for k = 1:N_g]
u⁻_lattices_by_g = [Array{Array{Float64, 3},1}(undef, M_amp) for k = 1:N_g]
u⁺_xy_lattices_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
u⁻_xy_lattices_by_g = [zeros(Float64, L₁, L₂) for k = 1:N_g]
u⁺_avg_by_g = [Array{Float64, 1}(undef, M) for k = 1:N_g]
u⁻_avg_by_g = [Array{Float64, 1}(undef, M) for k = 1:N_g]

# We preform M measurements
t_meas = @elapsed @time for m = 1:M
    dE_lists = nMCSEnUp!(cubs, Δt)[1]
    δE = [sum(dE_list) for dE_list in dE_lists]

    # Sort results in terms of temperature and extract ρˣˣₖ₂.
    for i = 1:N_g
        E_by_g[i][m] = E_prev[i] + δE[i]
        E_prev[i] += δE[i]
    end

    # Measure vortices and structure function
    for k = 1:N_g
        proj_V⁺, proj_V⁻ = xyVortexSnapshot(cubs[k])
        V⁺_avg_by_g[k] .+= proj_V⁺; V⁻_avg_by_g[k] .+= proj_V⁻
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        S⁺_by_g[k][m], S⁻_by_g[k][m] = structureFunction(proj_V⁺, proj_V⁻)

        proj_Vx, proj_Vy = xyVortexSnapshotXYBasis(cubs[k])
        Vx_avg_by_g[k] .+= proj_Vx; Vy_avg_by_g[k] .+= proj_Vy
        Sx, Sy = structureFunction(proj_Vx, proj_Vy)
        Sx_avg_by_g[k] .+= Sx; Sy_avg_by_g[k] .+= Sy
    end

    # Measure amplitudes
    for k = 1:N_g
        u⁺_lattice, u⁻_lattice = chiralAmplitudeSnapshot(cubs[k])
        # Save z-averaged lattices
        u⁺_avg_z = avgZ(u⁺_lattice); u⁻_avg_z = avgZ(u⁻_lattice)
        u⁺_xy_lattices_by_g[k] .+= u⁺_avg_z
        u⁻_xy_lattices_by_g[k] .+= u⁻_avg_z
        u⁺_avg_by_g[k][m] = mean(u⁺_avg_z)
        u⁻_avg_by_g[k][m] = mean(u⁻_avg_z)
        # For the first M_amp measurements we save the entire matrix as well
        if m <= M_amp
            u⁺_lattices_by_g[k][m] = u⁺_lattice
            u⁻_lattices_by_g[k][m] = u⁻_lattice
        end
    end

end
V⁺_avg_by_g = V⁺_avg_by_g./M; V⁻_avg_by_g  = V⁻_avg_by_g./M
Vx_avg_by_g = Vx_avg_by_g./M; Vy_avg_by_g = Vy_avg_by_g./M
Sx_avg_by_g = Sx_avg_by_g./M; Sy_avg_by_g = Sy_avg_by_g./M
final_lattices = [getLattice(cub) for cub in cubs]
vortices_by_g = [vortexSnapshot(lattice, getSyst(cubs[1])) for lattice in final_lattices]

println("Measurements used $(round(t_meas/3600; digits=1)) h, which means it took on average $(round(t_meas/(M*N_g); digits=1)) s pr. measurement.")
println("and $(round(t_meas/(M*Δt*N*N_g)*1e6; digits=2)) μs pr. MCS pr. lattice site.")



# Saving results to file
################################################################################################
for k = 1:N_g
    JLD.save(out_folders[k]*"/energies.jld", "Es", E_by_g[k])

    JLD.save(out_folders[k]*"/meta.jld", "L1", L₁, "L2", L₂, "L3", L₃, "M", M, "M_amp", M_amp, "dt", Δt, "T", T, "n", n, "m", m, "kappa", κ₅, "nu", ν, "g", gs[k])
    JLD.save(out_folders[k]*"/vorticity.jld", "vortexes", vortices_by_g[k], "sp", S⁺_by_g[k], "sm", S⁻_by_g[k])
    JLD.save(out_folders[k]*"/real_vorticity.jld", "vp_avg", V⁺_avg_by_g[k], "vm_avg", V⁻_avg_by_g[k])
    JLD.save(out_folders[k]*"/XY_vorticity.jld", "vx_avg", Vx_avg_by_g[k], "vy_avg", Vy_avg_by_g[k], "sx_avg", Sx_avg_by_g[k], "sy_avg", Sy_avg_by_g[k])
    JLD.save(out_folders[k]*"/amplitudes.jld", "up_lattices", u⁺_lattices_by_g[k], "um_lattices", u⁻_lattices_by_g[k], "up_xy", u⁺_xy_lattices_by_g[k], "um_xy", u⁻_xy_lattices_by_g[k], "up_avg", u⁺_avg_by_g, "um_avg", u⁻_avg_by_g)
end

JLD.save("final_state_nu=$(round(ν; digits=3))_L=$(L)_T=$(T).jld", "lattices", final_lattices, "T", T, "n", n, "m", m, "kappa", κ₅, "nu", ν, "gs", gs, "controls", [getControls(cub) for cub in cubs])
