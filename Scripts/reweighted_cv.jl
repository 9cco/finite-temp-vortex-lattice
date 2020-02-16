# This script takes as input JLD files in the calling directory that starts with "energy_table".
# It is assumed these files contain a dictionary with an "Ens_by_T" and "T_list" entry that
# contains an Array{Array{Float64, 1}} of energies and Array{Float64, 1} of corresponding temperatures.
#
# Using this data, the script calculates the average specific heat with jackknife error,
# assuming a system size of N and using num_blocks blocks.
#
# It then uses the data for temperatures between T_min and T_max to reweigh the specific heat to the
# temperatures set in rw_Ts. All raw and reweighted cvs are saved to JLD data file.
#
include("./functions_post.jl")

num_blocks = 2^4
T_min = 1.7
T_max = 1.85
N = 42^3
rw_Ts = vcat([1.85, 1.825, 1.8, 1.775, 1.75, 1.725, 1.7],
             [x for x in 1.7:0.01:1.85])


rw_Ts = sort(unique(rw_Ts))
rw_βs = [1/T for T in rw_Ts]

# Idiot check
if T_min > minimum(rw_Ts) || T_max < maximum(rw_Ts);  println("Warning: We are reweighting outside the range of input temperatures"); end

println("Reweighting specific heat cv using energy measurements for temperatures between $(T_min) - $(T_max),
for a system with N = $(N), using $(num_blocks) bins. Reweighting at $(length(rw_Ts)) new temperatures.")

raw_Es, raw_Ts = loadEnergyFiles()
N_T = length(raw_Ts)
raw_βs = [1/T for T in raw_Ts]

# Plotting all the raw data we have available
println("Data from $(N_T) temperatures found. Calculating Cv for all raw data.")
cv_raw_avg_by_T, cv_raw_err_by_T = calculateCvUncertainty(raw_Es, raw_βs, N; num_blocks=num_blocks)
raw_T_min = minimum(raw_Ts)
raw_T_max = maximum(raw_Ts)

res_Ts, cv_avgs, cv_errs, rw_cvs, rw_cv_errs = reweightAndErrorCvForInterval(rw_Ts, raw_Ts, raw_Es, (T_min, T_max), N)

JLD.save("rw_cv_data_$(num_blocks)blk.jld", "raw_Ts", raw_Ts, "raw_cv_avgs", cv_raw_avg_by_T, "raw_cv_errs", cv_raw_err_by_T, "rw_Ts", rw_Ts, "rw_cvs", rw_cvs, "rw_cv_errs", rw_cv_errs, "Ts", res_Ts, "cv_avgs", cv_avgs, "cv_errs", cv_errs, "num_blocks", num_blocks)


