function initiateStates(init_files::Array{String, 1}, L₁::I, L₂::I, L₃::I, g::P, ν::P, κ₅::P, κs::Array{P, 1}, n::I2, m::I2, 
                        split::Tuple{I3, I3, I3}; T_start=8.0) where {I<:Int, P<:Real, I2<:Int, I3<:Int}
    N_κ = length(κs); workers_pr_state = Int(prod(split))

    println("\nINITIATING STATES\n#############################################################")
    #
    # We check if the staging folder (could be current folder) includes an initial state file for
    # each state going to be simulated.
    has_init_files = true; for path in init_files; if !isfile(path); has_init_files=false; break; end; end;
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

    cubs, systs, init_temps
end

function estimateRuntime!(cubs::Array{Cuboid, 1}, N::I, M_est::I2, M_col::I2, M_th::I2, M::I2, Δt::I2) where {I<:Int, I2<:Int}
    N_c = length(cubs)

    # Do M_est MC-sweeps on first
    println("Estimating runtime")
    t_meas = @elapsed for i = 1:M_est
        nMCSEnUp!(cubs, Δt)

        En_res = [energy(cub) for cub in cubs]
        measureStates(cubs; full_amplitude_lattice=true)
    end
    t_meas = t_meas/M_est
    t_MCS = t_meas/(Δt*length(cubs))
    println("We did $(M_est) measurement steps using average $(round(t_meas; digits=1)) s pr measurement. 
    This yields $(round(t_MCS*1000; digits=1)) ms on average for each MCS")
    println("and $(round(t_MCS/N*1e6; digits=1)) μs on average pr. lattice site.")

    # Estimating ETC
    println("Measurements and thermalization will do $(Δt*M + M_th + M_col) MCS for each of the $(N_c)
    states, which has an ETC of $(round((Δt*M + M_th + M_col)*N_c*t_MCS/3600; digits=2)) h")
    print("Continue with thermalization and measurements? (y/n): ")
    flush(stdout)
    nothing
end

function cooldownStates!(cubs::Array{Cuboid, 1}, N::I, M_col::I2, N_steps::I2, M_est::I2, init_temps::Array{R, 1}, T::R,
                         out_folders::Array{String, 1}) where {I<:Int, I2<:Int, R<:Real}
    N_c = length(cubs)
    println("\nStarted cooldown at: $(Dates.format(now(), "HH:MM"))")

    # Do M_th MCS in N_steps steps
    M_pr_step = ceil(Int64, M_col/N_steps)
    M_col = M_pr_step*N_steps
    # Each step is associated with a different set of temperatures given by each row in the matrix
    temp_mt = genGeometricTemperatureSteps(init_temps, fill(T, N_c), N_steps)
    # Energy storage
    E_matrix = Array{Float64, 2}(undef, M_col+1, N_c);
    for (i, cub) = enumerate(cubs); E_matrix[1,i] = energy(cub); end
    # Acceptance rate storage
    accepts_matrix = Array{Int64, 2}(undef, M_col, N_c);

    println("Cooling down from T=$(init_temps) to T=$(T) using $(N_steps) steps with $(M_pr_step) MCS pr. step.")
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

    println("Cooldown used $(round(t_col/3600; digits=1)) h, i.e. $(round(t_col/(M_col*N_c*N)*1e6; digits=2)) μs pr. lattice site")
    E_matrix = E_matrix./N;

    # Saving Cooldown to separate ν-folders
    for k = 1:N_c
        JLD.save(out_folders[k]*"/cooldown.jld", "E_list", E_matrix[:, k], "accepts_list", accepts_matrix[:,k], "temp_list", temp_mt[:, k], "N_steps", N_steps, "M_pr_step", M_pr_step, "M_est", M_est)
    end

    M_col, t_col
end

function thermalizeStates!(cubs::Array{Cuboid, 1}, M_th::I, M_est::I, M_col::I, out_folders::Array{String, 1}) where I<:Int
    N_c = length(cubs)
    println("Thermalizing at target temperatures for additional $(M_th) MCS")
    println("\nStarted thermalization at: $(Dates.format(now(), "HH:MM"))")
    flush(stdout)

    E_therm = Array{Float64, 2}(undef, M_th, N_c);
    for j = 1:N_c; E_therm[1, j] = energy(cubs[j]); end
    t_th = @elapsed δE_lists = nMCSEnUp!(cubs, M_th-1)[1]
    for i = 2:M_th
        for k = 1:N_c
            E_therm[i, k] = E_therm[i-1, k] + δE_lists[k][i-1]
        end
    end

    E_therm = E_therm./N

    # Update time-estimation
    t_MCS1 = t_th/(M_th*N_c)

    therm_lattices = [getLattice(cub) for cub in cubs]

    # Saving thermalization
    for k = 1:N_c
        JLD.save(out_folders[k]*"/thermalization.jld", "e_thm", E_therm[:,k], "syst", systs[k], "lattice", therm_lattices[k], "M_th", M_th, "M_est", M_est, "M_col", M_col)
    end
    println("Thermalization used $(round(t_th/3600; digits=1)) h, i.e. $(round(t_MCS1/(N_c*N)*1e6; digits=2)) μs pr. lattice site. Energies and states saved to file.")

    t_MCS1
end

function measureStates(cubs::Array{Cuboid, 1}; full_amplitude_lattice=false)
    N_c = length(cubs)

    # Temp Storages
    proj_V⁺ = Array{Array{Float64, 2}, 1}(undef, N_c)
    proj_V⁻ = Array{Array{Float64, 2}, 1}(undef, N_c)
    S⁺ = Array{Array{Float64, 2}, 1}(undef, N_c)
    S⁻ = Array{Array{Float64, 2}, 1}(undef, N_c)
    proj_Vx = Array{Array{Float64, 2}, 1}(undef, N_c)
    proj_Vy = Array{Array{Float64, 2}, 1}(undef, N_c)
    Sx = Array{Array{Float64, 2}, 1}(undef, N_c)
    Sy = Array{Array{Float64, 2}, 1}(undef, N_c)
    if full_amplitude_lattice
        u⁺_lattices = Array{Array{Float64, 3}, 1}(undef, N_c)
        u⁻_lattices = Array{Array{Float64, 3}, 1}(undef, N_c)
    end
    u⁺_avg_z = Array{Array{Float64, 2}, 1}(undef, N_c)
    u⁻_avg_z = Array{Array{Float64, 2}, 1}(undef, N_c)
    u⁺_avg = Array{Float64, 1}(undef, N_c)
    u⁻_avg = Array{Float64, 1}(undef, N_c)
    δu² = Array{Float64, 1}(undef, N_c)

    for k = 1:N_c
        # Measure vortices and structure function
        proj_V⁺[k], proj_V⁻[k] = xyVortexSnapshot(cubs[k])
        # Find the Fourier transform of the projected vortices ∀ k ∈ k_matrix
        S⁺[k], S⁻[k] = structureFunction(proj_V⁺[k], proj_V⁻[k])
        
        proj_Vx[k], proj_Vy[k] = xyVortexSnapshotXYBasis(cubs[k])
        Sx[k], Sy[k] = structureFunction(proj_Vx[k], proj_Vy[k])
    
        # Measure amplitudes
        u⁺_lattice, u⁻_lattice = chiralAmplitudeSnapshot(cubs[k])
        u⁺_avg_z[k] = avgZ(u⁺_lattice); u⁻_avg_z[k] = avgZ(u⁻_lattice)
        u⁺_avg[k] = mean(u⁺_avg_z[k])
        u⁻_avg[k] = mean(u⁻_avg_z[k])
        if full_amplitude_lattice
            u⁺_lattices[k] = u⁺_lattice
            u⁻_lattices[k] = u⁻_lattice
        end
        δu²[k] = mean(u⁺_lattice.^2 .- u⁻_lattice.^2)
    end
    
    # Measure phase twist derivatives
    dH_01s = firstDerivativeTwist(cubs, 0, 1)
    dH_10s = firstDerivativeTwist(cubs, 1, 0)
    dH_11s = firstDerivativeTwist(cubs, 1, 1)
    d²H_01s = secondDerivativeTwist(cubs, 0, 1)
    d²H_10s = secondDerivativeTwist(cubs, 1, 0)
    d²H_11s = secondDerivativeTwist(cubs, 1, 1)

    if full_amplitude_lattice
        return (proj_V⁺, proj_V⁻, S⁺, S⁻, proj_Vx, proj_Vy, Sx, Sy, u⁺_lattices, u⁻_lattices, u⁺_avg_z, u⁻_avg_z, u⁺_avg, u⁻_avg, δu²,
            dH_01s, dH_10s, dH_11s, d²H_01s, d²H_10s, d²H_11s)
    else
        return (proj_V⁺, proj_V⁻, S⁺, S⁻, proj_Vx, proj_Vy, Sx, Sy, u⁺_avg_z, u⁻_avg_z, u⁺_avg, u⁻_avg, δu²,
            dH_01s, dH_10s, dH_11s, d²H_01s, d²H_10s, d²H_11s)
    end
end
