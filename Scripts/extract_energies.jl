using JLD

folders = filter(s -> isdir(s) && occursin(r"C4_model.*_mult_kap", s), readdir())
N_F = length(folders)
home_folder = pwd()

# Storage
Ens_by_T = Array{Array{Float64, 1}, 1}(undef, N_F)
T_list = Array{Float64, 1}(undef, N_F)

# Enter first folder and collect metadata
#cd(folders[1])

cd(home_folder)

for (i_f, folder) = enumerate(folders)
    cd(folder)
    println("Extracting energies in $(folder)")

    # Saving to storage

    ens_di = JLD.load("energies.jld")
    Es = ens_di["Es"]
    meta_di = JLD.load("meta.jld")
    T = meta_di["T"]

    T_list[i_f] = T
    Ens_by_T[i_f] = Es

    cd(home_folder)
end

T_min = minimum(T_list)
T_max = maximum(T_list)

# Saving table information in jld file

out_file = "energy_tableT=$(T_min)-$(T_max).jld"

if isfile(out_file)
    println("Warning: $(out_file) already exists. Overwrite? (y/n): ")
    user_input = readline(stdin)
    if occursin(r"[nN]", user_input)
        exit()
    end
end

JLD.save(out_file, "T_list", T_list, "Ens_by_T", Ens_by_T)
