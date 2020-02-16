include("./functions_post.jl")

plot_folder = "Temp_Series"

folders = filter(s -> isdir(s) && occursin(r"C4_model.*_lowField", s), readdir())
num_blocks = 2^6

#S⁺_points = Array{Tuple{Int64, Int64}, 1}(undef, 0)
N_F = length(folders)
home_folder = pwd()

# Enter first folder and collect metadata
cd(folders[1])
L₁, L₂, L₃, N, n, m, f, g, ν, κ₅ = extractMetaSyst()
cd(home_folder)

S⁺_points = [(Int(L₁/2)+3, 1+Int(L₂/2)-5) # Left hex
             ,(Int(L₁/2)+6, 1+Int(L₂/2)+1) # Upper right hex
             ,(Int(L₁/2)+3, 1+Int(L₂/2)+6) # Right Hex
             ,(Int(L₁/2)+5, 1+Int(L₂/2)-5) # Upper left quad
             ,(Int(L₁/2)+5, 1+Int(L₂/2)+5) # Upper right quad
             ,(Int(L₁/2)+5, 1+Int(L₂/2)) # Mid top ring
             ,(Int(L₁/2)+20, 1+Int(L₂/2)+20) # Plasma point
             ,(Int(L₁/2), 1+Int(L₂/2)) # Mid point
#             ,(Int(L₁/2)-5, 5+Int(L₂/2)) # Quad point
#             ,(Int(L₁/2)+5, 4+Int(L₂/2)) # Stripe point
]
#S⁺_points = Array{Tuple{Int64, Int64}, 1}(undef, L₁*L₂)
#for x = 1:L₁
#    for y = 1:L₂
#        S⁺_points[(x-1)*L₂+y] = (x, y)
#    end
#end
point_names = ["hex_ul", "hex_ur", "hex_right",
            "quad_ul", "quad_ur", "mid_top", "plasma", "mid"]
hex_lim = 0.2
quad_lim = 0.006
N_p = length(S⁺_points)


# Now we extract all measurements that we want to plot
Ts, Ess, S⁺s_Ts_points, u⁺ss, u⁻ss, δu²s = loadMeasurements(folders, S⁺_points; accu=false, fuzzy=false)
N_T = length(Ts)
println("Loaded $(length(u⁺ss)) amp measurements and total $(length(Ts)) temperatures.")
println("Calculate jackknife averages and errors")
mkcd(plot_folder)

# Calculate jackknife errors and averages.
E_avgs, E_errs = jackknifeMeanErr(Ess./N; num_blocks=num_blocks)
u⁺_avgs, u⁺_errs = jackknifeMeanErr(u⁺ss; num_blocks=num_blocks)
u⁻_avgs, u⁻_errs = jackknifeMeanErr(u⁻ss; num_blocks=num_blocks)
δu²_avgs, δu²_errs = jackknifeMeanErr(δu²s; num_blocks=num_blocks)
normalization = (two_pi*f*L₁*L₂)^2
S⁺_points_avgs = [Array{Float64, 1}(undef, N_T) for p = 1:N_p]
S⁺_points_errs = [Array{Float64, 1}(undef, N_T) for p = 1:N_p]
for p = 1:N_p
    for k = 1:N_T
        S⁺_points_avgs[p][k], S⁺_points_errs[p][k] = jackknifeMeanErr(S⁺s_Ts_points[p][k]./normalization; num_blocks=num_blocks)
    end
end
# For the first point, we plot the variance plot by bins
singleSeriesVariancePlot(S⁺s_Ts_points[1][1], num_blocks; filename="variance_point_$(getRelativeTuple(S⁺_points[1], L₁,L₂))_T=$(Ts[1])_$(num_blocks)blk.pdf")
println("Calculated $(N_p) structure function points")
cv_avgs, cv_errs = calculateCvUncertainty(Ess, [1/T for T in Ts], N; num_blocks=num_blocks)

T_min = minimum(Ts)
T_max = maximum(Ts)

##########################################################################################################################

# Plotting 
margin = 0.005
title="g = $(round(g; digits=2)), ν = $(round(ν; digits=2)), fL₁ = $(round(f*L₁; digits=2)), L₁=$(L₁)"

u⁺_avgs, u⁺_errs = jackknifeMeanErr(u⁺ss; num_blocks=num_blocks)
u⁻_avgs, u⁻_errs = jackknifeMeanErr(u⁻ss; num_blocks=num_blocks)
δu²_avgs, δu²_errs = jackknifeMeanErr(δu²s; num_blocks=num_blocks)
#[L"\langle\bar{u}^+\rangle", L"\langle\bar{u}^-\rangle"]
println("Plot temperature series")
pltTempSeries(Ts, cv_avgs, cv_errs; ylabel=L"\frac{C_v}{N}", title=title, filename="cv_by_T.pdf")
pltTempSeries(Ts, E_avgs, E_errs; ylabel="E", title=title, filename="E_by_T.pdf")
pltTempSeries(Ts[1:length(u⁺_avgs)], [u⁺_avgs, u⁻_avgs], [u⁺_errs, u⁻_errs]; ylabel=L"\bar{|u^\pm|}", labels=["u+", "u-"], title=title, filename="amplitudes_by_T.pdf")
pltTempSeries(Ts[1:length(δu²_avgs)], δu²_avgs, δu²_errs; ylabel=L"\langle\overline{u^{+2}-u^{-2}}\rangle", title=title, filename="du_by_T.pdf")
flux_lattice_folder = "Fourier"
if !isdir(flux_lattice_folder)
    mkdir(flux_lattice_folder)
end

plot_ending = "by_T.pdf"
# Plotting hex-points
for p = 1:3
    rel_point = getRelativeTuple(S⁺_points[p], L₁, L₂)
    pltTempSeries(Ts, S⁺_points_avgs[p], S⁺_points_errs[p]; ylabel=L"S^+(K)", ylims=(0, hex_lim), title="K = 2π[$(rel_point[1]), $(rel_point[2])]/L", filename=flux_lattice_folder*"/S+_$(point_names[p])_$(rel_point)_$(plot_ending)")
#    pltTempSeries(Ts, S⁺_avg_Ts_points[p]./normalization; ylabel=L"S^+(K)", title="K = 2π$(S⁺_points[p])/L ", filename="S+_avg_$(S⁺_points[p])_by_T.pdf")
end
# Plotting quad points
for p = 4:5
    rel_point = getRelativeTuple(S⁺_points[p], L₁, L₂)
    pltTempSeries(Ts, S⁺_points_avgs[p], S⁺_points_errs[p]; ylabel=L"S^+(K)", ylims=(0, quad_lim), title="K = 2π[$(rel_point[1]), $(rel_point[2])]/L", filename=flux_lattice_folder*"/S+_$(point_names[p])_$(rel_point)_$(plot_ending)")
#    pltTempSeries(Ts, S⁺_avg_Ts_points[p]./normalization; ylabel=L"S^+(K)", title="K = 2π$(S⁺_points[p])/L ", filename="S+_avg_$(S⁺_points[p])_by_T.pdf")
end
# Plotting rest
for p = 6:N_p
    rel_point = getRelativeTuple(S⁺_points[p], L₁, L₂)
#    if maximum(S⁺_points_avgs[p]) > 0.02
pltTempSeries(Ts, S⁺_points_avgs[p], S⁺_points_errs[p]; ylabel=L"S^+(K)", title="K = 2π[$(rel_point[1]), $(rel_point[2])]/L", filename=flux_lattice_folder*"/S+_$(rel_point)_$(plot_ending)")
#    end
#    pltTempSeries(Ts, S⁺_avg_Ts_points[p]./normalization; ylabel=L"S^+(K)", title="K = 2π$(S⁺_points[p])/L ", filename="S+_avg_$(S⁺_points[p])_by_T.pdf")
end



# Plotting energy histograms

plotEnergyHistograms(Ess, Ts)

# Saving table information in jld file

JLD.save("table_infoT=$(T_min)-$(T_max).jld", "Ts", Ts, "cv_avgs", cv_avgs, "cv_errs", cv_errs)


cd("../")
