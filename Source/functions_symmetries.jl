# functions_symmetries.jl
# Contains functionality to measure order parameters in the U(1) and chiral sectors



###################################################################################################
#       U(1) gauge stiffness used to measure the superconducting phase transition
#
###################################################################################################

#--------------------------------------------------------------------------------------------------
# Dual gauge stiffness for some small fourrier mode k (typically (0, 2π/l, 0). Implement xx, yy and
# zz components separately
#--------------------------------------------------------------------------------------------------
# zz-component gauge stiffness
function gaugeStiffnessZZ{T<:Real}(k::Array{T,1}, ψ::State)
    sumGS = Complex(0.0)
    L = ψ.consts.L
    if (L != ψ.consts.L₃)
        throw(error("The geometry is not cubic, L ≂̸ L₃"))
    end
    
    for h_pos=1:L, v_pos=1:L, z_pos=1:L
        r = [h_pos-1, L-v_pos, L-z_pos] #Spør 9cco om dette, tror det blir riktig
        ϕ = ψ.lattice[v_pos, h_pos, z_pos]
        ϕᵣ₊₁ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₁
        ϕᵣ₊₂ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₂
        sumGS += (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])*exp(im*(k⋅r))
    end
    return (abs(sumGS))^2 / (2π*L)^3
end

# xx-component gauge stiffness
function gaugeStiffnessXX{T<:Real}(k::Array{T,1}, ψ::State)
    sumGS = Complex(0.0)
    L = ψ.consts.L
    if (L != ψ.consts.L₃)
        throw(error("The geometry is not cubic, L ≂̸ L₃"))
    end
    
    for h_pos=1:L, v_pos=1:L, z_pos=1:L
        r = [h_pos-1, L-v_pos, L-z_pos] #Spør 9cco om dette, tror det blir riktig
        ϕ = ψ.lattice[v_pos, h_pos, z_pos]
        ϕᵣ₊₂ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₂
        ϕᵣ₊₃ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₃
        sumGS += (ϕ.A[2] + ϕᵣ₊₂.A[3] - ϕᵣ₊₃.A[2] - ϕ.A[3])*exp(im*(k⋅r))
    end
    return (abs(sumGS))^2 / (2π*L)^3
end

# yy-component gauge stiffness
function gaugeStiffnessYY{T<:Real}(k::Array{T,1}, ψ::State)
    sumGS = Complex(0.0)
    L = ψ.consts.L
    if (L != ψ.consts.L₃)
        throw(error("The geometry is not cubic, L ≂̸ L₃"))
    end
    
    for h_pos=1:L, v_pos=1:L, z_pos=1:L
        r = [h_pos-1, L-v_pos, L-z_pos] #Spør 9cco om dette, tror det blir riktig
        ϕ = ψ.lattice[v_pos, h_pos, z_pos]
        ϕᵣ₊₁ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₁
        ϕᵣ₊₃ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₃
        sumGS += (ϕ.A[3] + ϕᵣ₊₃.A[1] - ϕᵣ₊₁.A[3] - ϕ.A[1])*exp(im*(k⋅r))
    end
    return (abs(sumGS))^2 / (2π*L)^3
end
#--------------------------------------------------------------------------------------------------
# For a state, make M measurements of the gaugestiffness components, sampling every Δt MCS
function gaugeStiffnessMeasure!(ψ::State, sim::Controls, M::Int64, Δt::Int64)
    L = ψ.consts.L
    
    #Itialialise measurement arrays and fourrier modes
    ρˣˣₖ₂ = Array{Float64}(M)
    ρˣˣₖ₃ = Array{Float64}(M)
    ρʸʸₖ₃ = Array{Float64}(M)
    ρʸʸₖ₁ = Array{Float64}(M)
    ρᶻᶻₖ₁ = Array{Float64}(M)
    ρᶻᶻₖ₂ = Array{Float64}(M)
    k₁ = [2π/L, 0.0, 0.0]
    k₂ = [0.0, 2π/L, 0.0]
    k₃ = [0.0, 0.0, 2π/L]
    u⁺_array = Array{Float64}(M)
    u⁻_array = Array{Float64}(M)
    A_array = Array{Float64}(M)
    E_array = Array{Float64}(M)
    θ_array = Array{Float64}(M)
    θu_array = Array{Float64}(M)
    
    #Initial measurement
    ρˣˣₖ₂[1] = gaugeStiffnessXX(k₂, ψ)
    ρˣˣₖ₃[1] = gaugeStiffnessXX(k₃, ψ)
    ρʸʸₖ₃[1] = gaugeStiffnessYY(k₃, ψ)
    ρʸʸₖ₁[1] = gaugeStiffnessYY(k₁, ψ)
    ρᶻᶻₖ₁[1] = gaugeStiffnessZZ(k₁, ψ)
    ρᶻᶻₖ₂[1] = gaugeStiffnessZZ(k₂, ψ)
    E_array[1] = E(ψ)
    

    A_abs = 0.0 
    u⁺ = 0.0
    u⁻ = 0.0
    θ = 0.0
    θu = 0.0
    for x=1:L, y=1:L, z=1:L
        ϕ = ψ.lattice[x,y,z]
        ϕᵣ₊₁ = ψ.nb[x,y,z].ϕᵣ₊₁
        ϕᵣ₊₂ = ψ.nb[x,y,z].ϕᵣ₊₂
        ϕᵣ₊₃ = ψ.nb[x,y,z].ϕᵣ₊₃
        u⁺ += ψ.lattice[x,y,z].u⁺^2
        u⁻ += ψ.lattice[x,y,z].u⁻^2
        A_abs += ((ψ.lattice[x,y,z].A[1])^2 + (ψ.lattice[x,y,z].A[2])^2
                    + (ψ.lattice[x,y,z].A[3])^2 )
        θ += (cos(ϕ.θ⁺-ϕᵣ₊₁.θ⁺) + cos(ϕ.θ⁺-ϕᵣ₊₂.θ⁺) + cos(ϕ.θ⁺-ϕᵣ₊₃.θ⁺))
        θu += ϕ.u⁺*(ϕᵣ₊₁.u⁺*cos(ϕ.θ⁺-ϕᵣ₊₁.θ⁺) + ϕᵣ₊₂.u⁺*cos(ϕ.θ⁺-ϕᵣ₊₂.θ⁺) + ϕᵣ₊₃.u⁺*cos(ϕ.θ⁺-ϕᵣ₊₃.θ⁺))
    end
    u⁺_array[1] = u⁺/(L^3 )
    u⁻_array[1] = u⁻/(L^3 )
    A_array[1] = A_abs/(L^3 )
    θ_array[1] = θ/(3*L^3)
    θu_array[1] = θu/(3*L^3)

    #M-1 remaining measurements
    for m = 2:M
        #Make Δt MCS
        for mcs = 1:Δt
            mcSweep!(ψ, sim)
        end
         
        #Meaurements
        ρˣˣₖ₂[m] = gaugeStiffnessXX(k₂, ψ)
        ρˣˣₖ₃[m] = gaugeStiffnessXX(k₃, ψ)
        ρʸʸₖ₃[m] = gaugeStiffnessYY(k₃, ψ)
        ρʸʸₖ₁[m] = gaugeStiffnessYY(k₁, ψ)
        ρᶻᶻₖ₁[m] = gaugeStiffnessZZ(k₁, ψ)
        ρᶻᶻₖ₂[m] = gaugeStiffnessZZ(k₂, ψ)
        E_array[m] = E(ψ)
        
        A_abs = 0.0 
        u⁺ = 0.0
        u⁻ = 0.0
        θ = 0.0
        θu = 0.0
        for x=1:L, y=1:L, z=1:L
            ϕ = ψ.lattice[x,y,z]
            ϕᵣ₊₁ = ψ.nb[x,y,z].ϕᵣ₊₁
            ϕᵣ₊₂ = ψ.nb[x,y,z].ϕᵣ₊₂
            ϕᵣ₊₃ = ψ.nb[x,y,z].ϕᵣ₊₃
            u⁺ += ψ.lattice[x,y,z].u⁺
            u⁻ += ψ.lattice[x,y,z].u⁻
            A_abs += ((ψ.lattice[x,y,z].A[1])^2 + (ψ.lattice[x,y,z].A[2])^2
                    + (ψ.lattice[x,y,z].A[3])^2 )
            θ += (cos(ϕ.θ⁺-ϕᵣ₊₁.θ⁺) + cos(ϕ.θ⁺-ϕᵣ₊₂.θ⁺) + cos(ϕ.θ⁺-ϕᵣ₊₃.θ⁺))
            θu += ϕ.u⁺*(ϕᵣ₊₁.u⁺*cos(ϕ.θ⁺-ϕᵣ₊₁.θ⁺) + ϕᵣ₊₂.u⁺*cos(ϕ.θ⁺-ϕᵣ₊₂.θ⁺) + ϕᵣ₊₃.u⁺*cos(ϕ.θ⁺-ϕᵣ₊₃.θ⁺))
        end
        u⁺_array[m] = u⁺/(L^3 )
        u⁻_array[m] = u⁻/(L^3 )
        A_array[m] = A_abs/(L^3 )
        θ_array[m] = θ/(3*L^3)
        θu_array[m] = θu/(3*L^3)
    end
    return (ρˣˣₖ₂, ρˣˣₖ₃, ρʸʸₖ₃, ρʸʸₖ₁, ρᶻᶻₖ₁, ρᶻᶻₖ₂, u⁺_array, u⁻_array, A_array, E_array,
            θ_array, θu_array)
end
    
    

#--------------------------------------------------------------------------------------------------
# Take a list of uncorrelated, thermalised states and perform M measurements on these split over
# np workers and a master process. Measurements are taken every Δt steps, and only the averages are
# saved. In this first version we only measure the gauge stiffness.
function parallelMeasureGS!(ψ_list::Array{State,1}, sim::Controls, M::Int64, Δt::Int64)
    syst = ψ_list[1].consts
    L = syst.L
    
    #Initialise parallel process
    np = nprocs()
    length(ψ_list) >= np || throw(error("Not enough states in list"))
    M_min = Int(floor(M/np))    #Minimum amount of measurements per worker
    nw = M % np                 #Number of workers doing one extra measurement
    futures = [Future() for i=1:(np-1)]
    
    println("Starting $(M) measurements on $(np) processes on an $(L)x$(L)x$(syst.L₃) system,
            correspoding to $(M_min + Int(ceil(nw/np))) measurements and 
            $((M_min+ceil(Int64, nw/np))*Δt) MCS pr. process")
    
    #Start the +1 workers
    for i=1:nw
        futures[i] = @spawn gaugeStiffnessMeasure!(ψ_list[i], sim, M_min+1, Δt)
    end
    #Remaining workers
    for i=1:np-nw-1
        futures[nw+i] = @spawn gaugeStiffnessMeasure!(ψ_list[i], sim, M_min, Δt)
    end
    #Master process
    (ρˣˣₖ₂, ρˣˣₖ₃, ρʸʸₖ₃, ρʸʸₖ₁, ρᶻᶻₖ₁, ρᶻᶻₖ₂, u⁺, u⁻, A, E_arr,
        θ_arr, θu_arr) = gaugeStiffnessMeasure!(ψ_list[np], sim, M_min, Δt) 
   
    
    println("Measurements completed, collecting results from processes")
    #Gather results from workers
    for i=1:np-1
        (newρˣˣₖ₂, newρˣˣₖ₃, newρʸʸₖ₃, newρʸʸₖ₁, newρᶻᶻₖ₁, newρᶻᶻₖ₂, newu⁺, newu⁻, newA, newE,
            newθ, newθu) = fetch(futures[i])
        ρˣˣₖ₂ = vcat(ρˣˣₖ₂, newρˣˣₖ₂)
        ρˣˣₖ₃ = vcat(ρˣˣₖ₃, newρˣˣₖ₃)
        ρʸʸₖ₃ = vcat(ρʸʸₖ₃, newρʸʸₖ₃)
        ρʸʸₖ₁ = vcat(ρʸʸₖ₁, newρʸʸₖ₁)
        ρᶻᶻₖ₁ = vcat(ρᶻᶻₖ₁, newρᶻᶻₖ₁)
        ρᶻᶻₖ₂ = vcat(ρᶻᶻₖ₂, newρᶻᶻₖ₂)
        u⁺ = vcat(u⁺, newu⁺)
        u⁻ = vcat(u⁻, newu⁻)
        A = vcat(A,newA)
        E_arr = vcat(E_arr, newE)
        θ_arr = vcat(θ_arr, newθ)
        θu_arr = vcat(θu_arr, newθu)
    end
    
    println("Processing results")
    #Calculate averages
    AVρˣˣₖ₂ = mean(ρˣˣₖ₂)
    AVρˣˣₖ₃ = mean(ρˣˣₖ₃)
    AVρʸʸₖ₃ = mean(ρʸʸₖ₃)
    AVρʸʸₖ₁ = mean(ρʸʸₖ₁)
    AVρᶻᶻₖ₁ = mean(ρᶻᶻₖ₁)
    AVρᶻᶻₖ₂ = mean(ρᶻᶻₖ₂)
    AVu⁺ = mean(u⁺)
    AVu⁻ = mean(u⁻)
    AVA = mean(A)
    AVE = mean(E_arr)/L^3
    AVθ = (mean(θ_arr))^2
    AVθu = (mean(θu_arr))^2
    
    
    #Calculate the standard error 
    #TODO: Here we might also have to account for the autocorr time
    SEρˣˣₖ₂ = std(ρˣˣₖ₂)/sqrt(M)
    SEρˣˣₖ₃ = std(ρˣˣₖ₃)/sqrt(M)
    SEρʸʸₖ₃ = std(ρʸʸₖ₃)/sqrt(M)
    SEρʸʸₖ₁ = std(ρʸʸₖ₁)/sqrt(M)
    SEρᶻᶻₖ₁ = std(ρᶻᶻₖ₁)/sqrt(M)
    SEρᶻᶻₖ₂ = std(ρᶻᶻₖ₂)/sqrt(M)
    
    return (AVρˣˣₖ₂, SEρˣˣₖ₂, AVρˣˣₖ₃, SEρˣˣₖ₃, AVρʸʸₖ₃, SEρʸʸₖ₃, AVρʸʸₖ₁, SEρʸʸₖ₁,
            AVρᶻᶻₖ₁, SEρᶻᶻₖ₁, AVρᶻᶻₖ₂, SEρᶻᶻₖ₂, AVu⁺, AVu⁻, AVA, u⁺, AVE, AVθ, AVθu)
end
#--------------------------------------------------------------------------------------------------   
#Functionality to measure the helicity modulus
function helicityModulusContributions(ψ::State)
    sumCosX = 0.0
    sumCosY = 0.0
    sumCosZ = 0.0
    sumSinX = 0.0
    sumSinY = 0.0
    sumSinZ = 0.0

    L = ψ.consts.L
    for h_pos = 1:L, v_pos = 1:L, z_pos=1:L
        ϕ = ψ.lattice[h_pos, v_pos, z_pos]
        ϕᵣ₊₁ = ψ.nb[h_pos, v_pos, z_pos].ϕᵣ₊₁
        ϕᵣ₊₂ = ψ.nb[h_pos, v_pos, z_pos].ϕᵣ₊₂
        ϕᵣ₊₃ = ψ.nb[h_pos, v_pos, z_pos].ϕᵣ₊₃
        sumCosX += ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺)
        sumCosY += ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺)
        sumCosZ += ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺ - ϕ.θ⁺)
        sumSinX += ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺)
        sumSinY += ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺)
        sumSinZ += ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺ - ϕ.θ⁺)

    end
    return (sumCosX, sumCosY, sumCosZ, sumSinX^2, sumSinY^2, sumSinZ^2)
end


