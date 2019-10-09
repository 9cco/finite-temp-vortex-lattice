# Current intention:
# Figure out if we can use 1/4 of thermalization steps and still get a valid state, while using new code
# with measured helicity modulus and variable MGT term. We would also like, when this is imported over
# to Fram to run two L=64 states using a single script, where each state is divided in 16 subcubes,
# instead of a single state divided on 32 subcubes as previously.

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

@everywhere src_path = "/home/nicolai/Documents/Work/PhD/Numerikk/MC/finite-temp-vortex-lattice/Source/Grid/"
#@everywhere src_path = "/cluster/home/fredrkro/mc/Source/Grid/"
out_path = "/home/nicolai/mc/Scripts/"
#out_path = "/cluster/home/fredrkro/mc/Data/ExperimentalRuns/"
staging_path = ""
#staging_path = "/cluster/home/fredrkro/Staging/"
@everywhere push!(LOAD_PATH, $src_path)
@everywhere using CuboidModule

@everywhere include(src_path*"observables.jl")
@everywhere include(src_path*"utilities.jl")


using JLD


# Enter data directory for structure function
fixRC()
# We run a simulation with parameter κ supplied by script
if length(ARGS) != 1
    println("ERROR: Need κ supplied as script argument.")
    exit(-1)
end
ν = parse(Float64, ARGS[1])             # Fermi surface anisotropy
κs = [1.0]                         # MGT weights
g = 0.316227766                               # Gauge couplings
κ₅ = 1.0
N_κ = length(κs)

# Other parameters
M_est = 2^8#2^7  # Number of MC-steps to do at T_start before cooldown.
M_col = 2^17 # Number of MC-steps to use for cooling down the systems
N_steps = 2^10   # Number of temperatures to go through from high temp before reaching the final temperature. (will divide M_col) must be >= 2
M_th = 2^17  # Number of MC-steps to do for thermalization.
M = 2^12     # Number of measurements
M_amp = 2^5   # Number of these measurements that will be amplitude measurements, i.e. we need M >= M_amp
Δt = 2^5      # Number of MC-steps between each measurement
# L is assumed to be even.
L = 32     # System length
L₁ = L
L₂ = L
L₃ = L
split = (1,1,4)#(1,1,2)#
N = L₁*L₂*L₃
# Specify extended landau gauge using n,m s.t. constant Aᵢ(Lⱼ+1) = Aᵢ(0) mod 2π
# n | L₁ and m | L₂
n = 1; m = 0
f = n/L₁ - m/L₂
T = 2.0
T_start = 2*4.0+0.1   # Start-temperature of states when thermalizing. 


# Print parameters
println("\nStarted script at: $(Dates.format(now(), "HH:MM"))")
println("Target temp: $(T)")
println("f = $(f)")
println("n = $(n)")
println("m = $(m)")
println("L = $(L)")
println("ν = $(ν)")
println("κ₅ = $(κ₅)")
println("κs = $(κs)")
println("g = $(g)")
if f < 0
    println("Number of anti-vortices: $(-f*L₁*L₂)")
else
    println("Number of vortices: $(f*L₁*L₂)")
end
flush(stdout)
cd(out_path)
out_folders = ["full_model_L=$(L)_T=$(round(T; digits=2))_g=$(round(g; digits=3))_nu=$(round(ν; digits=3))_kap=$(κ)_n=$(n)_m=$(m)_mult_kap" for κ in κs]

# Checking available workers
workers_pr_state = Int(split[1]*split[2]*split[3])
needed_workers = workers_pr_state*N_κ + N_κ - 1
if nprocs()-1 < needed_workers
    println("ERROR: Only $(nprocs()-1) workers have been allocated. $(needed_workers) are needed to split all
of the $((needed_workers-N_κ+1)/workers_pr_state) states into $(split[1]*split[2]*split[3]) subcuboids and then
run these in parallel using $(N_κ) additional workers.")
    exit()
end



# Initiating states
################################################################################################
println("\nINITIATING STATES\n#############################################################")
#
# We check if the staging folder (could be current folder) includes an initial state file for
# each state going to be simulated.
init_files = [staging_path*"init_state_g=$(round(g; digits=3))_nu=$(round(ν; digits=3))_kap=$(round(κ, digits=3)).jld" for κ in κs]
has_init_files = true; for path in init_files; if !isfile(path); global has_init_files=false; break; end; end;
if has_init_files

    # Construct initial states
    cubs = Array{Cuboid, 1}(undef, N_κ)
    systs = Array{SystConstants, 1}(undef, N_κ)
    init_temps = Array{Float64, 1}(undef, N_κ)

    for i = 1:N_κ
        # Making control system to compare with
        control_syst = SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,κs[i],n,m)
        cubs[i], systs[i], init_temps[i] = constructStateFromJLD(init_files[i], control_syst, split, (i-1)*workers_pr_state + N_κ + 1)
    end
    
    println("$(N_κ) state(s) created from initiation file with κ-value(s):")
    println(κs)
    println("Controls of read state:")
    printControls(cubs)
    
    # Tuning initial controls
    flush(stdout)
    AR_est = [estimateAR!(cub)[1] for cub in cubs]
    println("Current AR with these controls:")
    for (k, est) = enumerate(AR_est); println("AR(g=$(round(κs[k]; digits=1))): ≈ $(round(est*100; digits=2))%"); end

else # If no initial file is found, then we contruct initial states at temperature T_start and thermalize them
    # Make ab inito un-correlated phases state if no initial_state is given.
    systs = [SystConstants(L₁,L₂,L₃,1/g^2,ν,κ₅,κ,n,m) for κ in κs]
    cubs = Array{Cuboid, 1}(undef, N_κ)
    for k = 1:N_κ
        cubs[k] = Cuboid(6, systs[k], split, 1/T_start; u⁺=1/√(2), u⁻=1/√(2), pid_start=(k-1)*workers_pr_state + N_κ + 1)
    end
    init_temps = [T_start for syst in systs]

    println("No initialization file found. Thermalizing ab-inito states.")
    flush(stdout)
    # Tune and thermalize start-temperature
    t_init_th = @elapsed nMCS!(cubs, ceil(Int64, M_th/4))

    println("Adjusting simulation constants in initial states.")
    for k = 1:N_κ; setUpdates!(cubs[k]; θ_max = 3.13, u_max = 1.0/√(2), A_max = 1e-1); end # If we set updates manually
    println("Controls after adjustment:")
    printControls(cubs)
    flush(stdout)

    t_init_th2 = @elapsed nMCS!(cubs, ceil(Int64, 3*M_th/4))
    println("Thermalization of ab-inito states used $(round((t_init_th+t_init_th2)/3600; digits=1)) h i.e.
$(round((t_init_th+t_init_th2)/(M_th*N)*1e6; digits=1)) μs on average pr. lattice site.")
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

    En_res = [energy(cub) for cub in cubs]
    for cub in cubs; chiralAmplitudeSnapshot(cub); xyVortexSnapshot(cub); end
    
    firstDerivativeTwist(cubs, 0, 1)
    firstDerivativeTwist(cubs, 1, 0)
    firstDerivativeTwist(cubs, 1, 1)
    secondDerivativeTwist(cubs, 0, 1)
    secondDerivativeTwist(cubs, 1, 0)
    secondDerivativeTwist(cubs, 1, 1)
end
t_meas = t_meas/M_est
t_MCS = t_meas/(Δt*length(cubs))
println("We did $(M_est) measurement steps using average $(round(t_meas; digits=1)) s pr measurement. 
This yields $(round(t_MCS*1000; digits=1)) ms on average for each MCS")
println("and $(round(t_MCS/N*1e6; digits=1)) μs on average pr. lattice site.")

# Estimating ETC
println("Measurements and thermalization will do $(Δt*M + M_th + M_col) MCS for each of the $(N_κ)
g-values, which has an ETC of $(round((Δt*M + M_th + M_col)*N_κ*t_MCS/3600; digits=2)) h")
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
temp_mt = genGeometricTemperatureSteps(init_temps, fill(T, N_κ), N_steps)
# Energy storage
E_matrix = Array{Float64, 2}(undef, M_col+1, N_κ);
for (i, cub) = enumerate(cubs); E_matrix[1,i] = energy(cub); end
E_therm = Array{Float64, 2}(undef, M_th, N_κ);
# Acceptance rate storage
accepts_matrix = Array{Int64, 2}(undef, M_col, N_κ);

println("Cooling down from T=$(init_temps) using $(N_steps) steps with $(M_pr_step) MCS pr. step.")
flush(stdout)

t_col = @elapsed for step = 1:N_steps
    # First set the temperatures associated with this step to each system.
    for (k, cub) = enumerate(cubs); setTemp!(cub, temp_mt[step, k]); end

    # Do the MCS associated with this step
    δE_lists, accepts_lists = nMCSEnUp!(cubs, M_pr_step)
    for i = 1:M_pr_step
        for (k, cub) = enumerate(cubs)
            δE = δE_lists[k][i]
            accepts = accepts_lists[k][i]
            E_matrix[(step-1)*M_pr_step+i+1, k] = E_matrix[(step-1)*M_pr_step+i, k] + δE
            accepts_matrix[(step-1)*M_pr_step+i, k] = accepts
        end
    end

    # In the end of a step we correct energy rounding-error
    for (j, cub) = enumerate(cubs); E_matrix[step*M_pr_step+1, j] = energy(cub); end
end

println("Cooldown used $(round(t_col/3600; digits=1)) h, i.e. $(round(t_col/(M_col*N_κ*N)*1e6; digits=2)) μs pr. lattice site")
E_matrix = E_matrix./N;

# Saving Cooldown to separate ν-folders
for k = 1:N_κ
    JLD.save(out_folders[k]*"/cooldown.jld", "E_list", E_matrix[:, k], "accepts_list", accepts_matrix[:,k], "temp_list", temp_mt[:, k], "N_steps", N_steps, "M_pr_step", M_pr_step, "M_est", M_est)
end

t_MCS0 = t_col/(M_col*N_κ)



# Thermalization
################################################################################################
println("Thermalizing at target temperatures for additional $(M_th) MCS")
println("\nStarted thermalization at: $(Dates.format(now(), "HH:MM"))")
flush(stdout)

for j = 1:N_κ; E_therm[1, j] = energy(cubs[j]); end
t_th = @elapsed δE_lists = nMCSEnUp!(cubs, M_th-1)[1]
for i = 2:M_th
    for k = 1:N_κ
        E_therm[i, k] = E_therm[i-1, k] + δE_lists[k][i-1]
    end
end
AR_est = [estimateAR!(cub)[1] for cub in cubs]

E_therm = E_therm./N

# Update time-estimation
t_MCS1 = t_th/(M_th*N_κ)
t_MCS_av = (t_MCS0 + t_MCS1)/2

therm_lattices = [getLattice(cub) for cub in cubs]

# Saving thermalization
for k = 1:N_κ
    JLD.save(out_folders[k]*"/thermalization.jld", "e_thm", E_therm[:,k], "syst", systs[k], "lattice", therm_lattices[k], "M_th", M_th, "M_est", M_est, "M_col", M_col)
end
println("Thermalization used $(round(t_th/3600; digits=1)) h, i.e. $(round(t_MCS1/(N_κ*N)*1e6; digits=2)) μs pr. lattice site. Energies and states saved to file.")
println("After thermalization the estimated AR are:")
for (k, est) = enumerate(AR_est); println("AR(g=$(round(κs[k]; digits=1))): ≈ $(round(est*100; digits=2))%"); end



# Doing measurements
################################################################################################
println("\nMEASUREMENTS\n#############################################################")
println("We will now do $(Δt*M) MCS pr ν value  which gives
$(M) measurements and has an ETC of $(round(Δt*M*N_κ*t_MCS_av/3600; digits=2)) h")
println("Started measurements at: $(Dates.format(now(), "HH:MM"))")
flush(stdout)

# Setup storage for dual stifness and energy
E_by_κ = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
E_prev = [energy(cub) for cub in cubs]
# Storage for vortices and dual vortices
S⁺_by_κ = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_κ]
S⁻_by_κ = [Array{Array{Float64, 2},1}(undef, M) for k = 1:N_κ]
V⁺_avg_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
V⁻_avg_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
Vx_avg_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
Vy_avg_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
Sx_avg_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
Sy_avg_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
# Storage for amplitude matrices
u⁺_lattices_by_κ = [Array{Array{Float64, 3},1}(undef, M_amp) for k = 1:N_κ]
u⁻_lattices_by_κ = [Array{Array{Float64, 3},1}(undef, M_amp) for k = 1:N_κ]
u⁺_xy_lattices_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
u⁻_xy_lattices_by_κ = [zeros(Float64, L₁, L₂) for k = 1:N_κ]
u⁺_avg_by_κ = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
u⁻_avg_by_κ = [Array{Float64, 1}(undef, M) for k = 1:N_κ]

# Storage for helicity modulus derivatives
dH_01_x = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_01_y = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_01_z = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_10_x = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_10_y = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_10_z = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_11_x = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_11_y = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
dH_11_z = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_01_x = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_01_y = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_01_z = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_10_x = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_10_y = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_10_z = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_11_x = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_11_y = [Array{Float64, 1}(undef, M) for k = 1:N_κ]
d²H_11_z = [Array{Float64, 1}(undef, M) for k = 1:N_κ]


# We preform M measurements
t_meas = @elapsed @time for m = 1:M
    dE_lists = nMCSEnUp!(cubs, Δt)[1]
    δE = [sum(dE_list) for dE_list in dE_lists]

    # Sort results in terms of temperature.
    for i = 1:N_κ
        E_by_κ[i][m] = E_prev[i] + δE[i]
        E_prev[i] += δE[i]
    end

    # Measure vortices and structure function
    for k = 1:N_κ
        proj_V⁺, proj_V⁻ = xyVortexSnapshot(cubs[k])
        V⁺_avg_by_κ[k] .+= proj_V⁺; V⁻_avg_by_κ[k] .+= proj_V⁻
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        S⁺_by_κ[k][m], S⁻_by_κ[k][m] = structureFunction(proj_V⁺, proj_V⁻)

        proj_Vx, proj_Vy = xyVortexSnapshotXYBasis(cubs[k])
        Vx_avg_by_κ[k] .+= proj_Vx; Vy_avg_by_κ[k] .+= proj_Vy
        Sx, Sy = structureFunction(proj_Vx, proj_Vy)
        Sx_avg_by_κ[k] .+= Sx; Sy_avg_by_κ[k] .+= Sy
    end

    # Measure amplitudes
    for k = 1:N_κ
        u⁺_lattice, u⁻_lattice = chiralAmplitudeSnapshot(cubs[k])
        # Save z-averaged lattices
        u⁺_avg_z = avgZ(u⁺_lattice); u⁻_avg_z = avgZ(u⁻_lattice)
        u⁺_xy_lattices_by_κ[k] .+= u⁺_avg_z
        u⁻_xy_lattices_by_κ[k] .+= u⁻_avg_z
        u⁺_avg_by_κ[k][m] = mean(u⁺_avg_z)
        u⁻_avg_by_κ[k][m] = mean(u⁻_avg_z)
        # For the first M_amp measurements we save the entire matrix as well
        if m <= M_amp
            u⁺_lattices_by_κ[k][m] = u⁺_lattice
            u⁻_lattices_by_κ[k][m] = u⁻_lattice
        end
    end

    # Measure phase twist derivatives
    dH_01s = firstDerivativeTwist(cubs, 0, 1)
    dH_10s = firstDerivativeTwist(cubs, 1, 0)
    dH_11s = firstDerivativeTwist(cubs, 1, 1)
    d²H_01s = secondDerivativeTwist(cubs, 0, 1)
    d²H_10s = secondDerivativeTwist(cubs, 1, 0)
    d²H_11s = secondDerivativeTwist(cubs, 1, 1)
    for k = 1:N_κ
        dH_01_x[k][m], dH_01_y[k][m], dH_01_z[k][m] = dH_01s[k]
        dH_10_x[k][m], dH_10_y[k][m], dH_10_z[k][m] = dH_10s[k] 
        dH_11_x[k][m], dH_11_y[k][m], dH_11_z[k][m] = dH_11s[k]

        d²H_01_x[k][m], d²H_01_y[k][m], d²H_01_z[k][m] = d²H_01s[k] 
        d²H_10_x[k][m], d²H_10_y[k][m], d²H_10_z[k][m] = d²H_10s[k]
        d²H_11_x[k][m], d²H_11_y[k][m], d²H_11_z[k][m] = d²H_11s[k]
    end


end
V⁺_avg_by_κ = V⁺_avg_by_κ./M; V⁻_avg_by_κ  = V⁻_avg_by_κ./M
Vx_avg_by_κ = Vx_avg_by_κ./M; Vy_avg_by_κ = Vy_avg_by_κ./M
Sx_avg_by_κ = Sx_avg_by_κ./M; Sy_avg_by_κ = Sy_avg_by_κ./M
final_lattices = [getLattice(cub) for cub in cubs]
vortices_by_κ = [vortexSnapshot(lattice, getSyst(cubs[1])) for lattice in final_lattices]

println("Measurements used $(round(t_meas/3600; digits=1)) h, which means it took on average $(round(t_meas/(M*N_κ); digits=1)) s pr. measurement.")
println("and $(round(t_meas/(M*Δt*N*N_κ)*1e6; digits=2)) μs pr. MCS pr. lattice site.")



# Saving results to file
################################################################################################
for k = 1:N_κ
    JLD.save(out_folders[k]*"/energies.jld", "Es", E_by_κ[k])

    JLD.save(out_folders[k]*"/meta.jld", "L1", L₁, "L2", L₂, "L3", L₃, "M", M, "M_amp", M_amp, "dt", Δt, "T", T, "n", n, "m", m, "kap5", κ₅, "kap", κs[k], "nu", ν, "g", g)
    JLD.save(out_folders[k]*"/vorticity.jld", "vortexes", vortices_by_κ[k], "sp", S⁺_by_κ[k], "sm", S⁻_by_κ[k])
    JLD.save(out_folders[k]*"/real_vorticity.jld", "vp_avg", V⁺_avg_by_κ[k], "vm_avg", V⁻_avg_by_κ[k])
    JLD.save(out_folders[k]*"/XY_vorticity.jld", "vx_avg", Vx_avg_by_κ[k], "vy_avg", Vy_avg_by_κ[k], "sx_avg", Sx_avg_by_κ[k], "sy_avg", Sy_avg_by_κ[k])
    JLD.save(out_folders[k]*"/amplitudes.jld", "up_lattices", u⁺_lattices_by_κ[k], "um_lattices", u⁻_lattices_by_κ[k], "up_xy", u⁺_xy_lattices_by_κ[k], "um_xy", u⁻_xy_lattices_by_κ[k], "up_avg", u⁺_avg_by_κ, "um_avg", u⁻_avg_by_κ)

    # Saving Helicity moduli
    JLD.save(out_folders[k]*"/hel_mod.jld", "dH_01_x", dH_01_x[k], "dH_01_y", dH_01_y[k], "dH_01_z", dH_01_z[k],
             "dH_10_x", dH_10_x[k], "dH_10_y", dH_10_y[k], "dH_10_z", dH_10_z[k],
             "dH_11_x", dH_11_x[k], "dH_11_y", dH_11_y[k], "dH_11_z", dH_11_z[k],
             "d2H_01_x", d²H_01_x[k], "d2H_01_y", d²H_01_y[k], "d2H_01_z", d²H_01_z[k],
             "d2H_10_x", d²H_10_x[k], "d2H_10_y", d²H_10_y[k], "d2H_10_z", d²H_10_z[k],
             "d2H_11_x", d²H_11_x[k], "d2H_11_y", d²H_11_y[k], "d2H_11_z", d²H_11_z[k])

    # Saving final state
    JLD.save(out_folders[k]*"/final_state_g=$(round(g; digits=3))_nu=$(round(ν; digits=3))_kap=$(round(κs[k], digits=3)).jld", "lattice", final_lattices[k], "syst", getSyst(cubs[k]), "T", getTemp(cubs[k]), "controls", getControls(cubs[k]))
end
