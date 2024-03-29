using Distributions
using Test
using StatsBase
using Plots
gr()


# Plots energy histograms for all the energy series contained in E_by_T and saves the plot in energy_histogram.pdf in the
# current folder unless another filename is specified. It uses automatic determination of bins for each energy-series
# unless nbins is specified.
# uses StatsBase, Plots - pyplot
# -----------------------------------------------------------------------------------------------------------
function plotEnergyHistograms(E_by_T::Vector{Vector{T}}, temps::Vector{T}; 
        nbins=-1, filename="energy_histograms.pdf") where {T}
    
    N₀ = length(E_by_T)
    title = "Energy histograms for $(length(temps)) temperatures"
    # If we are using automatic number of bins
    if nbins==-1
        # Make first histogram
        hist = fit(Histogram, E_by_T[1])
        edges = collect(hist.edges[1])
        xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
        ys = hist.weights
        plt = plot(xs, ys; label="T = $(round(temps[1]; digits=2))", xaxis="Energy pr. site", yaxis="number in bin", title=title);
        for k = 2:N₀
            hist = fit(Histogram, E_by_T[k])
            edges = collect(hist.edges[1])
            xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
            ys = hist.weights
            plot!(plt, xs, ys; label="T = $(round(temps[k]; digits=2))");
        end
    else
        # If we are using fixed number of bins
        # Make first histogram
        hist = fit(Histogram, E_by_T[1]; nbins=nbins)
        edges = collect(hist.edges[1])
        xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
        ys = hist.weights
        plt = plot(xs, ys; label="T = $(round(temps[1]; digits=2))", xaxis="Energy pr. site", yaxis="number in bin", title=title)
        for k = 2:N₀
            hist = fit(Histogram, E_by_T[k]; nbins=nbins)
            edges = collect(hist.edges[1])
            xs = [(edges[i] + edges[i+1])/2 for i = 1:length(edges)-1]
            ys = hist.weights
            plot!(plt, xs, ys; label="T = $(round(temps[k]; digits=2))");
        end
    end
    savefig(plt, filename);
    nothing
end

# -----------------------------------------------------------------------------------------------------------
# Plotting result of structureFunctionVortexLatticeAvg!
# -----------------------------------------------------------------------------------------------------------
# Plotting result of structureFunctionVortexLatticeAvg!
function plotStructureFunctionVortexLattice(avV⁺::Array{T, 2}, avV⁻::Array{T, 2}, 
        V⁺::Array{R, 2}, V⁻::Array{R,2}, avS⁺::Array{T, 2}, avS⁻::Array{T,2}, ks::Array{Array{T,1}, 2}) where {T<:Real, R<:Real}
    L = size(avV⁺,1)
    # Plotting structure factor
    L_k = size(ks, 1)
	k_x = [ks[L_k,x][1] for x=1:L_k]
    k_y = [ks[y,1][2] for y=1:L_k]
    plt = heatmap(k_x, k_y, avS⁺, title="average S+", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    display(plt)
    plt = heatmap(k_x, k_y, avS⁻, title="average S-", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    display(plt)
    
    # Removing middle-point
    println("S⁺(0) ≈ $(avS⁺[Int(ceil(L/2)), Int(ceil(1+L/2))])")
    temp⁺ = avS⁺[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))]
    avS⁺[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = avS⁺[Int(ceil(1+L_k/2)), Int(ceil(L_k/2))]
    println("S⁻(0) ≈ $(avS⁻[Int(ceil(L/2)), Int(ceil(1+L/2))])")
    temp⁻ = avS⁻[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))]
    avS⁻[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = avS⁻[Int(ceil(1+L_k/2)), Int(ceil(L_k/2))]
    
    # And then re-plotting
    plt = heatmap(k_x, k_y, avS⁺, 
        title="average S+ with S(0) removed", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    display(plt)
    plt = heatmap(k_x, k_y, avS⁻, 
        title="average S- with S(0) removed", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    display(plt)
    
    # Restoring middle point
    avS⁺[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = temp⁺
    avS⁻[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = temp⁻
    
    # Plotting vortex snapshots
    plt = heatmap(1:L, 1:L, V⁺, title="Snapshot of + component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    display(plt)
    plt = heatmap(1:L, 1:L, V⁻, title="Snapshot of - component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    display(plt)
    
    # Combining matrices
    combined_lattice = combineVortexLattices(V⁺, V⁻)
    plt = heatmap(1:L, 1:L, combined_lattice, title="Combination of snapshots", xlabel="x", ylabel="y", aspect_ratio=1)
    display(plt)
    
    # Finding the proportion of the different kinds of vortices in combined matrix
    av_vortex_kinds = zeros(9)
    for v_pos = 1:L, h_pos = 1:L
        if combined_lattice[v_pos, h_pos] >= 0
            av_vortex_kinds[combined_lattice[v_pos, h_pos]+1] += 1
        end
    end
    av_vortex_kinds *= 100/L^2
    println("The proportion of vortices (n⁺, n⁻) in snapshot")
    println("% of vortex kind (-1, -1): \t$(Int(round(av_vortex_kinds[1],0)))")
    println("% of vortex kind (-1, 0): \t$(Int(round(av_vortex_kinds[2],0)))")
    println("% of vortex kind (-1, 1): \t$(Int(round(av_vortex_kinds[3],0)))")
    println("% of vortex kind (0, -1): \t$(Int(round(av_vortex_kinds[4],0)))")
    println("% of vortex kind (0, 0): \t$(Int(round(av_vortex_kinds[5],0)))")
    println("% of vortex kind (0, 1): \t$(Int(round(av_vortex_kinds[6],0)))")
    println("% of vortex kind (1, -1): \t$(Int(round(av_vortex_kinds[7],0)))")
    println("% of vortex kind (1, 0): \t$(Int(round(av_vortex_kinds[8],0)))")
    println("% of vortex kind (1, 1): \t$(Int(round(av_vortex_kinds[9],0)))\n")
    
    # Calculating the sum of vortices of snapshot
    sum⁺ = 0.0
    sum⁻ = 0.0
    for h_pos = 1:L, v_pos = 1:L
        sum⁺ += V⁺[v_pos, h_pos]
        sum⁻ += V⁻[v_pos, h_pos]
    end
    println("Sum of + component vorticity in the snapshot: $(sum⁺)")
    println("Sum of - component vorticity in the snapshot: $(sum⁻)")
    flush(STDOUT)
    
    # Plotting vortex averages
    plt = heatmap(1:L, 1:L, avV⁺, title="Average + component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    display(plt)
    plt = heatmap(1:L, 1:L, avV⁻, title="Average - component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    display(plt)
end

# -----------------------------------------------------------------------------------------------------------
# Plotting and saving result of structureFunctionVortexLatticeAvg!
# Should be run inside of designated directory.
function plotStructureFunctionVortexLatticeS(ψ::State, avV⁺::Array{T, 2}, avV⁻::Array{T, 2}, 
        avS⁺::Array{T, 2}, avS⁻::Array{T,2}, ks::Array{Array{T,1}, 2}) where T<:Real
    L = size(avV⁺,1)
    # Plotting structure factor
    L_k = size(ks, 1)
	k_x = [ks[L_k,x][1] for x=1:L_k]
    k_y = [ks[y,1][2] for y=1:L_k]
    plt = heatmap(k_x, k_y, avS⁺, title="average S+", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    savefig(plt, "sfvl_avg_S+_plot.pdf")
    plt = heatmap(k_x, k_y, avS⁻, title="average S-", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    savefig(plt, "sfvl_avg_S-_plot.pdf")
    
    # Removing middle-point
    println("S⁺(0) ≈ $(avS⁺[Int(ceil(L/2)), Int(ceil(1+L/2))])")
    temp⁺ = avS⁺[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))]
    avS⁺[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = avS⁺[Int(ceil(1+L_k/2)), Int(ceil(L_k/2))]
    println("S⁻(0) ≈ $(avS⁻[Int(ceil(L/2)), Int(ceil(1+L/2))])")
    temp⁻ = avS⁻[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))]
    avS⁻[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = avS⁻[Int(ceil(1+L_k/2)), Int(ceil(L_k/2))]
    
    # And then re-plotting
    plt = heatmap(k_x, k_y, avS⁺, 
        title="average S+ with S(0) removed", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    savefig(plt, "sfvl_avg_S+_removed_plot.pdf")
    plt = heatmap(k_x, k_y, avS⁻, 
        title="average S- with S(0) removed", xlabel="k_x", ylabel="k_y", aspect_ratio=1)
    savefig(plt, "sfvl_avg_S-_removed_plot.pdf")
    
    # Restoring middle point
    avS⁺[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = temp⁺
    avS⁻[Int(ceil(L_k/2)), Int(ceil(1+L_k/2))] = temp⁻
    
    # Plotting vortex snapshots
    V⁺_mat, V⁻_mat = vortexSnapshot(ψ)
    L₃ = ψ.consts.L₃
    V⁺ = V⁺_mat[:,:,rand(1:L₃)]/two_pi; V⁻ = V⁻_mat[:,:,rand(1:L₃)]/two_pi
    plt = heatmap(1:L, 1:L, V⁺, title="Snapshot of + component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    savefig(plt, "sfvl_V+_snapshot_plot.pdf")
    plt = heatmap(1:L, 1:L, V⁻, title="Snapshot of - component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    savefig(plt, "sfvl_V-_snapshot_plot.pdf")
    
    # Combining matrices
    combined_lattice = combineVortexLattices(V⁺, V⁻)
    plt = heatmap(1:L, 1:L, combined_lattice, title="Combination of snapshots", xlabel="x", ylabel="y", aspect_ratio=1)
    savefig(plt, "sfvl_comb_snapshort_plot.pdf")
    
    # Finding the proportion of the different kinds of vortices in combined matrix
    av_vortex_kinds = zeros(9)
    for v_pos = 1:L, h_pos = 1:L
        if combined_lattice[v_pos, h_pos] >= 0
            av_vortex_kinds[combined_lattice[v_pos, h_pos]+1] += 1
        end
    end
    av_vortex_kinds *= 100/L^2
    println("The proportion of vortices (n⁺, n⁻) in snapshot")
    println("% of vortex kind (-1, -1): \t$(Int(round(av_vortex_kinds[1],0)))")
    println("% of vortex kind (-1, 0): \t$(Int(round(av_vortex_kinds[2],0)))")
    println("% of vortex kind (-1, 1): \t$(Int(round(av_vortex_kinds[3],0)))")
    println("% of vortex kind (0, -1): \t$(Int(round(av_vortex_kinds[4],0)))")
    println("% of vortex kind (0, 0): \t$(Int(round(av_vortex_kinds[5],0)))")
    println("% of vortex kind (0, 1): \t$(Int(round(av_vortex_kinds[6],0)))")
    println("% of vortex kind (1, -1): \t$(Int(round(av_vortex_kinds[7],0)))")
    println("% of vortex kind (1, 0): \t$(Int(round(av_vortex_kinds[8],0)))")
    println("% of vortex kind (1, 1): \t$(Int(round(av_vortex_kinds[9],0)))\n")
    
    # Calculating the sum of vortices of snapshot
    sum⁺ = 0.0
    sum⁻ = 0.0
    for h_pos = 1:L, v_pos = 1:L
        sum⁺ += V⁺[v_pos, h_pos]
        sum⁻ += V⁻[v_pos, h_pos]
    end
    println("Sum of + component vorticity in the snapshot: $(sum⁺)")
    println("Sum of - component vorticity in the snapshot: $(sum⁻)")
    flush(STDOUT)
    
    # Plotting vortex averages
    plt = heatmap(1:L, 1:L, avV⁺, title="Average + component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    savefig(plt, "sfvl_avg_V+_plot.pdf")
    plt = heatmap(1:L, 1:L, avV⁻, title="Average - component vorticity", xlabel="x", ylabel="y", aspect_ratio=1)
    savefig(plt, "sfvl_avg_V-_plot.pdf")
    
    # Saving average matrices to file
    writedlm("avg_S+.data", avS⁺, ":")
    writedlm("avg_S-.data", avS⁻, ":")
    writedlm("avg_V+.data", avV⁺, ":")
    writedlm("avg_V-.data", avV⁻, ":")
    writedlm("k_x_range.data", k_x, ":")
    writedlm("k_y_range.data", k_y, ":")
    
    return 1
end 


####################################################################################################################
#                            Diagnostic functions
#
####################################################################################################################


# -----------------------------------------------------------------------------------------------------------
# Plots various observables and fields of a state.
function diagnose(ψ::State; z=-1)
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    N = length(ψ.lattice)
    println("Checking state: $(checkState(ψ))")
    println("Checking neighbors: $(@test checkNeighbors(ψ))")
    println("Constants: g = $(sqrt(1/ψ.consts.g⁻²)), T = $(1/ψ.consts.β), L = $(ψ.consts.L), L₃ = $(ψ.consts.L₃)")
    println("Average energy: $(E(ψ)/N)")
    u⁺, u⁻ = meanAmplitudes(ψ)
    println("⟨u⁺⟩: $(u⁺)")
    Υˣ, Υʸ, Υᶻ = helicityModulus(ψ)
    println("Helicity modulus: $((Υˣ + Υʸ)/2)")
    println("⟨η⁺⟩: $(averageOrdParam(ψ)[1])")
    V⁺, V⁻ = vortexSnapshot(ψ)
    println("Sum of nᶻᵣ: $(sum(V⁺))")
    ρˣˣₖ₂, ρˣˣₖ₃, ρʸʸₖ₁, ρʸʸₖ₃, ρᶻᶻₖ₁, ρᶻᶻₖ₂ = gaugeStiffness([ψ])
    println("Dual stiffnesses:")
    println("ρˣˣ(k₂) =\t$(ρˣˣₖ₂[1])")
    println("ρˣˣ(k₃) =\t$(ρˣˣₖ₃[1])")
    println("ρʸʸ(k₁) =\t$(ρʸʸₖ₁[1])")
    println("ρʸʸ(k₃) =\t$(ρʸʸₖ₃[1])")
    println("ρᶻᶻ(k₁) =\t$(ρᶻᶻₖ₁[1])")
    println("ρᶻᶻ(k₂) =\t$(ρᶻᶻₖ₂[1])")
    

    # Find a random plane if not a valid plane has been specified.
    if z <= 0 || z > L₃
        z = rand(1:L₃)
    end
    println("\nLooking at plane $(z)")
    flush(STDOUT)
    phases = [ψ.lattice[v, h, z].θ⁺ for v = 1:L, h = 1:L]
    plt1 = heatmap(phases; aspect_ratio=1.0, title="θ⁺", xlabel="x", ylabel="y",
        colorbar_title="θ⁺", colorbar=:right, clims=(0,2π))
    amp = [ψ.lattice[v, h, z].u⁺ for v = 1:L, h = 1:L]
    plt2 = heatmap(amp; aspect_ratio=1.0, title="u⁺", xlabel="x", ylabel="y",
        colorbar_title="u⁺", colorbar=:right, clims=(0,1))
    plt = plot(plt1, plt2, layout=(1,2), size=(1000, 400))
    display(plt)
    Ax = [ψ.lattice[v, h, z].A[1] for v = 1:L, h = 1:L]
    Ay = [ψ.lattice[v, h, z].A[2] for v = 1:L, h = 1:L]
    Az = [ψ.lattice[v, h, z].A[3] for v = 1:L, h = 1:L]
    plt_x = heatmap(Ax; aspect_ratio=1.0, title="Aˣ", xlabel="x", ylabel="y", colorbar=:right)
    plt_y = heatmap(Ay; aspect_ratio=1.0, title="Aʸ", xlabel="x", ylabel="y",
        colorbar=:right)
    plt_z = heatmap(Az; aspect_ratio=1.0, title="Aᶻ", xlabel="x", ylabel="y",
        colorbar=:right)
    plt = plot(plt_x, plt_y, plt_z, layout=grid(1, 3, widths=[0.33, 0.33, 0.33]), legend=false, size=(1000,400))
    display(plt)
    flux_density_z = [fluxDensity(ψ, (v,h,z))[3] for v = 1:L, h = 1:L]
    plt1 = heatmap(flux_density_z; aspect_ratio=1.0, title="(∇×A)ᶻ", xlabel="x", ylabel="y", colorbar=:right)
    V⁺, V⁻ = vortexSnapshot(ψ)
    plt2 = heatmap(V⁺[:,:,z]./two_pi; aspect_ratio=1.0, title="nᶻᵣ", xlabel="x", ylabel="y", colorbar=:right, clims=(-1,1))
    plt = plot(plt1, plt2, layout=(1,2), size=(1000, 400))
    display(plt)
end

# -----------------------------------------------------------------------------------------------------------
# Testing how well our algorithms work for this choice of parameters.
function testSystem(syst::SystConstants, sim::Controls, M::Int64=2000)
    # Let's first look at the time it takes to equilibrate the state
    ψ = State(2, syst)
    ψ_old = copy(ψ)
    
    println("What happens after $(M) number of MCS")
    E_test = zeros(M)

    # Do mcSweep! M times
    for i = 1:M
        mcSweep!(ψ, sim)
        E_test[i] = E(ψ)
    end
    println("Checking if mcSweeped state has lower energy")
    println(@test E(ψ) < E(ψ_old))
    plt = plot(1:M, E_test, title="Development of internal free energy with Monte-Carlo sweeps", 
        xlabel="MCS", ylabel="E()")
    display(plt)
    
    # Then we find the proposal fraction for the update
    (res, stdev, fracs) = mcProposalFraction(ψ_old, M)
    plt = plot(1:M, fracs, title="Probability of update over time", ylabel="#accepted proposals / L^2", xlabel="MCS");
    display(plt)
    println("metropolisHastingUpdate! proposal was accepted $(round(res*100,1))±$(round(stdev*100,1))% of the times")
    
    # At infinite energy, metropolis Hasting will always accept the new state, therefore we expect
    # that after a couple of mcSweeps, the state will be completely changed
    println("Testing if mcSweep! gives completely different state when temperature is infinite, where completely different
    means that all values on all lattice sites are different. Thus we have proved that mcSweep! visits all lattice sites.")
    ψ = State(2, SystConstants(syst.L, syst.g⁻², syst.ν, syst.f, 0))
    ψ_old = copy(ψ)
    for i = 1:2
        mcSweep!(ψ)
    end
    println(@test latticeCompletelyDifferent(ψ,ψ_old))

    println("\nEquilibrium calculation\n----------------------------------------------------------------")
    ψ₁ = State(1, syst)
    ψ₂ = State(2, syst)
    println("Checking that a random state has lower energy than a completely correlated state")
    println(@test E(ψ₁) < E(ψ₂))
    flush(STDOUT)
    (t₀, dE, ψ₁, ψ₂) = findEquilibrium(syst, M)
    @show t₀
    T = size(dE,1)
    plt = plot(1:T,dE, title="Difference in energy for a correlated state - random state", xlabel="MCS", ylabel="E1-E2");
    display(plt)
end


####################################################################################################################
#                            Utility Functions
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Creates a new directory based on the simulation constants and enters it for further writing of files.
function mkcdSystemDirectory(syst::SystConstants, M::Int64, Δt::Int64)
    DATA_DIR = "Data"
    REMOTE_TREE = "/work/fredrkro/finite-temp-vortex-lattice/Data"
    LOCAL_TREE = "../Data"
    FOLDER_FILE = "last_folder.txt"
    acc = 3
    
    # Now we try to find the appropriate Data directory to put the new folder in.
    # Our first priority is checking if we are on the remote system.
    if ispath(REMOTE_TREE)
        println("It seems we are on a remote machine. Changing to appropriate directoy.")
        cd(REMOTE_TREE)
        # If the function is run from the Notebooks folder then we might have a data-directory as a neighbor
    elseif ispath(LOCAL_TREE)
        println("There is a neighboring Data directory. Changing to this.")
        cd(LOCAL_TREE)
        # If the current directory has a data directory below it we change to that.
    elseif isdir("./$(DATA_DIR)")
        println("A Data directory was found below the current. Changing.")
        cd(DATA_DIR)
    end
    
    # Make directory name based on system constants
    DIR_NAME = "3D_L_$(syst.L)_M_$(M)_T_$(round(1/syst.β,acc))_g_$(round(√(1/syst.g⁻²),acc))"
    DIR_NAME = DIR_NAME * "_NU_$(round(syst.ν,acc))_K_$(round(syst.κ₅,acc))_f_$(round(syst.f,acc))_T_$(round(1/syst.β,acc))"
    
    if isdir("./$(DIR_NAME)")
        print("The directory $(DIR_NAME) already exists.\nShould we delete it? [y/n]: ")
        answ = readline(STDIN)
        if answ == "n"
            return 0
        end
        rm("./$(DIR_NAME)", recursive=true)
		print("\n")
    end
    
    mkdir(DIR_NAME)
    println("Made directory at $(pwd())/$(DIR_NAME)")
    println("Writing folder name $(DIR_NAME) to $(FOLDER_FILE)")
    open(FOLDER_FILE, "w") do f
        write(f, "$(DIR_NAME)\n")
    end
    
    cd(DIR_NAME)
    return 1
end

# -----------------------------------------------------------------------------------------------------------
# In the current directory, writes the current system and simulations constants to a file system_values.data
function writeSimulationConstants(syst::SystConstants, sim::Controls, M::Int64, t₀::Int64, Δt::Int64, 
        filename::AbstractString = "system_values.data", dir::AbstractString="")
    if (dir != "")
        println("Making dir")
        makeDirRec(dir)
        dir_name = dir*"/"
    else
        dir_name = ""
    end
    open(dir_name*filename, "w") do f
        write(f, "L $(syst.L)\n")
        write(f, "Lz $(syst.L₃)\n")
        write(f, "g $(√(1/syst.g⁻²))\n")
        write(f, "NU $(syst.ν)\n")
        write(f, "K $(syst.κ₅)\n")
        write(f, "f $(syst.f)\n")
        write(f, "TEMP $(1/syst.β)\n")
        write(f, "INV_TEMP $(syst.β)\n")
        write(f, "NR_MEASUREMENTS $(M)\n")
        write(f, "MEASUREMENT_INTERVAL $(Δt)\n")
        write(f, "SIM_THETA_MAX $(sim.θmax)\n")
        write(f, "SIM_UMAX $(sim.umax)\n")
        write(f, "SIM_AMAX $(sim.Amax)\n")
        write(f, "THERM_T $(t₀)\n")
    end
end


# -----------------------------------------------------------------------------------------------------------
# Performing common initialization checks and plots for making sure the states are correctly
# thermalized.
# Saves to .pdf files in current directory.
function initializeParallelStatesS(syst::SystConstants, sim::Controls, dir::AbstractString="")
	# Parameters
	PLOT_FRACTION = 0.71
    PLOT_DE = 1e3

	# Parallel thermalization
	START_T = 1000 		# Number of MCS to start with
	EX = 1.8			# Factor to extend by if no dE<0 is found
	STD_FACTOR = 0.5	# How many standard errors the average have to be within.

	# Create initial states, reference at random and worker list from corrolation
    n_workers = max(nprocs()-1,1)
	ψ_ref = State(2, syst)
	ψ_w = [State(1,syst) for i = 1:n_workers]

	# Thermalize these states.
	println("Thermalizing $(n_workers+1) states")
	@time t₀, tᵢ, Tᵢ, E_ref, E_w, ψ_ref, ψ_w, sim_ref, sim_w = parallelThermalization!(ψ_ref, ψ_w, syst, sim, START_T, EX, STD_FACTOR)
    flush(STDOUT)
    if t₀ == -1
        throw(error("ERROR: Could not thermalize states."))
    end

	N = length(E_ref)
    dE_array = [E_ref - E_w[w,:] for w = 1:n_workers]
    # Plot from the point where the ref. energy and first worker is within PLOT_DE
    tₛ = 1
    for i = 1:N
        if abs(dE_array[1][i]) < PLOT_FRACTION
            tₛ = i
            break
        end
    end

	# Plot the last PLOT_FRACTION fraction of the energy intervals
    #tₛ = min(ceil(Int64, EX*START_T), N - floor(Int64, N*PLOT_FRACTION))
	int = tₛ:N
	n_workers == size(E_w,1) || throw(error("ERROR: Somehow the number of workers changed during thermalization"))

	# Converting the misc. energies into an array of arrays of energies
	energies = Array{Array{Float64, 1}, 1}(undef, n_workers+1)
	labels = Array{String, 1}(undef, n_workers+1)
	energies[n_workers+1] = E_ref[int]
	labels[n_workers+1] = "reference"
	for w = 1:n_workers
		energies[w] = E_w[w,int]
		labels[w] = "worker $(w)"
	end
    labels = reshape(labels, (1, n_workers+1))

    #Create directory to save plots. If dir is empty or "" it is saved in pwd
    if (dir != "")
        makeDirRec(dir)
        dir_name = dir*"/"
    else
        dir_name = ""
    end

	# Plotting the energies of reference and all workers
	plt = plot(int, energies, label=labels, xlabel="MCS", ylabel="Energy", title="Thermalization energies");
    savefig(plt, dir_name*"ini_equi_E_plot.pdf")
    
	# Plot energy difference for first worker
    plt = plot(int, dE_array[1][int], title="Energy difference", xlabel="MCS");
    savefig(plt, dir_name*"ini_equi_dE_plot.pdf")

	# Plot energy difference for all workers in averaging interval
    avg_int = tᵢ:Tᵢ
    dE_array = [dE_array[w][avg_int] for w = 1:n_workers]
	plt = plot(avg_int, dE_array, title="dE, average interval", xlabel="MCS", label=labels[1:n_workers], ylabel="En. diff.")
	savefig(plt, dir_name*"ini_equi_dE_avgint_plot.pdf")
    
	ψ₁ = ψ_w[1]
	ψ₂ = ψ_ref
	sim₁ = sim_w[1]
	sim₂ = sim_ref
    println("Calculating energies and acceptance rates")
    flush(STDOUT)
    T = 1000
    dE = zeros(2T)
    E₁ = zeros(2T)
    E₂ = zeros(2T)
    p₁ = zeros(2T)
    p₂ = zeros(2T)
    for i = 1:T
        E₁[i] = E(ψ₁)
        E₂[i] = E(ψ₂)
        dE[i] = E₂[i]-E₁[i]
        p₁[i] = mcSweepFrac!(ψ₁, sim₁)
        p₂[i] = mcSweepFrac!(ψ₂, sim₂)
    end
    adjustSimConstants!(sim₁, ψ₁)
    adjustSimConstants!(sim₂, ψ₂)
    for i = T+1:2T
        E₁[i] = E(ψ₁)
        E₂[i] = E(ψ₂)
        dE[i] = E₂[i]-E₁[i]
        p₁[i] = mcSweepFrac!(ψ₁, sim₁)
        p₂[i] = mcSweepFrac!(ψ₂, sim₂)
    end

	# Making %
	p₁ = 100 .* p₁
	p₂ = 100 .* p₂

	println("Saving plots to files")

    # Plotting results
	x_mcs = (t₀+1):(t₀+2T)
    plt = plot(x_mcs, dE, title="Energy difference", xlabel="MCS");
    savefig(plt, dir_name*"ini_extra_dE_plot.pdf")
    plt = plot(x_mcs, [E₁,E₂], title="Internal energies", xlabel="MCS");
    savefig(plt, dir_name*"ini_extra_E_plot.pdf")
    plt = plot(x_mcs, [p₁,p₂], title="Accept probabilities", xlabel="MCS", ylabel="%");
    savefig(plt, dir_name*"ini_extra_AR_plot.pdf")
    
	return (t₀, ψ_ref, sim_ref, ψ_w, sim_w)
end

# -----------------------------------------------------------------------------------------------------------
# Same as above, but now assuming only one process.
function initializeTwoStatesS(syst::SystConstants, sim::Controls)
    (t₀, E₁, E₂, dE, ψ₁, ψ₂, sim₁, sim₂) = findEquilibrium(syst, sim);
    flush(STDOUT)
    
    N = size(dE, 1)
    int = 1:N
    plt = plot(int, dE[int], title="Energy difference", xlabel="MCS");
    savefig(plt, "ini_equi_dE_plot.pdf")
    plt = plot(int, E₁[int], title="Internal energy 1", xlabel="MCS");
    savefig(plt, "ini_equi_E1_plot.pdf")
    plt = plot(int, E₂[int], title="Internal energy 2", xlabel="MCS");
    savefig(plt, "ini_equi_E2_plot.pdf")
    
    println("Performing extra MCS")
	last_pr = -1
    for i = 1:t₀
        mcSweep!(ψ₂)
        mcSweep!(ψ₁)
		this_pr = Int(round(i/t₀*100,0))
		# Only print to STDOUT when a new whole percentage is reached so that STDOUT is not
		# flooded and only print at 10% increments
		if (this_pr > last_pr && this_pr % 10 == 0)
			println("$(this_pr)%")
			last_pr = this_pr
		end
    end
    
    println("Calculating energies and acceptance rates")
    flush(STDOUT)
    T = 4000
    dE = zeros(2T)
    E₁ = zeros(2T)
    E₂ = zeros(2T)
    p₁ = zeros(2T)
    p₂ = zeros(2T)
    for i = 1:T
        E₁[i] = E(ψ₁)
        E₂[i] = E(ψ₂)
        dE[i] = E₂[i]-E₁[i]
        p₁[i] = mcSweepFrac!(ψ₁, sim₁)
        p₂[i] = mcSweepFrac!(ψ₂, sim₂)
    end
    adjustSimConstants!(sim₁, ψ₁)
    adjustSimConstants!(sim₂, ψ₂)
    for i = T+1:2T
        E₁[i] = E(ψ₁)
        E₂[i] = E(ψ₂)
        dE[i] = E₂[i]-E₁[i]
        p₁[i] = mcSweepFrac!(ψ₁, sim₁)
        p₂[i] = mcSweepFrac!(ψ₂, sim₂)
    end

	# Making %
	p₁ = 100 .* p₁
	p₂ = 100 .* p₂

    # Plotting results
	x_mcs = (2t₀+1):(2t₀+2T)
    plt = plot(1:2T, dE, title="Energy difference", xlabel="MCS", xticks=x_mcs);
    savefig(plt, "ini_extra_dE_plot.pdf")
    plt = plot(1:2T, E₁, title="Internal energy 1", xlabel="MCS", xticks=x_mcs);
    savefig(plt, "ini_extra_E1_plot.pdf")
    plt = plot(1:2T, E₂, title="Internal energy 2", xlabel="MCS", xticks=x_mcs);
    savefig(plt, "ini_extra_E2_plot.pdf")
    plt = plot(1:2T, p₁, title="Accept probability 1", xlabel="MCS", ylabel="%", xticks=x_mcs);
    savefig(plt, "ini_extra_AR1_plot.pdf")
    plt = plot(1:2T, p₂, title="Accept probability 2", xlabel="MCS", ylabel="%", xticks=x_mcs);
    savefig(plt, "ini_extra_AR2_plot.pdf")
    
    return (ψ₁, sim₁, ψ₂, sim₂, t₀)
end

#---------------------------------------------------------------------------------------------------
#Check if a directory exist, and make it if it does not
function makeDirRec(dir_name::AbstractString)
    folders = split(dir_name, "/")
    temp_dir = folders[1]
    if !isdir(temp_dir)
        mkdir(temp_dir)
    end
    
    for i=2:length(folders)
        temp_dir = temp_dir*"/"*folders[i]
        if !isdir(temp_dir)
            mkdir(temp_dir)
        end
    end
end
