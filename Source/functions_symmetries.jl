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
function gaugeStiffnessZZ(k::Array{T,1}, ψ::State) where T<:Real
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
        sumGS += (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])*exp(im*dot(k,r))
    end
    return (abs(sumGS))^2 / (two_pi^2*L^3)
end

# xx-component gauge stiffness
function gaugeStiffnessXX(k::Array{T,1}, ψ::State) where T<:Real
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
        sumGS += (ϕ.A[2] + ϕᵣ₊₂.A[3] - ϕᵣ₊₃.A[2] - ϕ.A[3])*exp(im*dot(k,r))
    end
    return (abs(sumGS))^2 / (two_pi^2*L^3)
end

# yy-component gauge stiffness
function gaugeStiffnessYY(k::Array{T,1}, ψ::State) where T<:Real
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
        sumGS += (ϕ.A[3] + ϕᵣ₊₃.A[1] - ϕᵣ₊₁.A[3] - ϕ.A[1])*exp(im*dot(k,r))
    end
    return (abs(sumGS))^2 / (two_pi^2*L^3)
end

# -------------------------------------------------------------------------------------------------
# Calculate gauge stiffness more effectively as a triplet
function gaugeStiffness(k::Array{T,1}, ψ::State) where T<:Real
    gs_xx = Complex(0.0); gs_yy = Complex(0.0); gs_zz = Complex(0.0)
    L = ψ.consts.L
    if (L != ψ.consts.L₃)
        throw(error("The geometry is not cubic, L ≂̸ L₃"))
    end

    for v_pos = 1:L, h_pos = 1:L, z_pos = 1:L
        r = [h_pos-1, L-v_pos, L-z_pos] #Spør 9cco om dette, tror det blir riktig
        ϕ = ψ.lattice[v_pos, h_pos, z_pos]
        ϕᵣ₊₁ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₁
        ϕᵣ₊₂ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₂
        ϕᵣ₊₃ = ψ.nb[v_pos, h_pos, z_pos].ϕᵣ₊₃
        gs_xx += (ϕ.A[2] + ϕᵣ₊₂.A[3] - ϕᵣ₊₃.A[2] - ϕ.A[3])*exp(im*dot(k,r))
        gs_yy += (ϕ.A[3] + ϕᵣ₊₃.A[1] - ϕᵣ₊₁.A[3] - ϕ.A[1])*exp(im*dot(k,r))
        gs_zz += (ϕ.A[1] + ϕᵣ₊₁.A[2] - ϕᵣ₊₂.A[1] - ϕ.A[2])*exp(im*dot(k,r))
    end

    norm = (2*π)^2*L^3
    return abs2(gs_xx)/norm, abs2(gs_yy)/norm, abs2(gs_zz)/norm
end

# -------------------------------------------------------------------------------------------------
# Given a list of states, calculate the gauge stiffnesses for each direction and each k, for each
# state in the list. Then return lists of the corresponding stiffnesses.
function gaugeStiffness(ψ_list::Array{State,1})
    M = length(ψ_list)
    L = ψ_list[1].consts.L

    ρˣˣₖ₂ = Array{Float64}(undef, M)
    ρˣˣₖ₃ = Array{Float64}(undef, M)
    ρʸʸₖ₁ = Array{Float64}(undef, M)
    ρʸʸₖ₃ = Array{Float64}(undef, M)
    ρᶻᶻₖ₁ = Array{Float64}(undef, M)
    ρᶻᶻₖ₂ = Array{Float64}(undef, M)
    k₁ = [2π/L, 0.0, 0.0]
    k₂ = [0.0, 2π/L, 0.0]
    k₃ = [0.0, 0.0, 2π/L]

    for m = 1:M
        ρˣˣₖ₂[m], _, ρᶻᶻₖ₂[m] = gaugeStiffness(k₂, ψ_list[m])
        ρˣˣₖ₃[m], ρʸʸₖ₃[m], _ = gaugeStiffness(k₃, ψ_list[m])
        _, ρʸʸₖ₁[m], ρᶻᶻₖ₁[m] = gaugeStiffness(k₁, ψ_list[m])
    end

    return ρˣˣₖ₂, ρˣˣₖ₃, ρʸʸₖ₁, ρʸʸₖ₃, ρᶻᶻₖ₁, ρᶻᶻₖ₂
end


#--------------------------------------------------------------------------------------------------
# For a state, make M measurements of the gaugestiffness components, sampling every Δt MCS
function gaugeStiffnessMeasure!(ψ::State, sim::Controls, M::Int64, Δt::Int64)
    L = ψ.consts.L
    
    #Itialialise measurement arrays and fourrier modes
    ρˣˣₖ₂ = Array{Float64}(undef, M)
    ρˣˣₖ₃ = Array{Float64}(undef, M)
    ρʸʸₖ₃ = Array{Float64}(undef, M)
    ρʸʸₖ₁ = Array{Float64}(undef, M)
    ρᶻᶻₖ₁ = Array{Float64}(undef, M)
    ρᶻᶻₖ₂ = Array{Float64}(undef, M)
    k₁ = [2π/L, 0.0, 0.0]
    k₂ = [0.0, 2π/L, 0.0]
    k₃ = [0.0, 0.0, 2π/L]
    u⁺_array = Array{Float64}(undef, M)
    u⁻_array = Array{Float64}(undef, M)
    A_array = Array{Float64}(undef, M)
    E_array = Array{Float64}(undef, M)
    ΥcosX = Array{Float64}(undef, M)
    ΥsinX = Array{Float64}(undef, M)
    ΥcosY = Array{Float64}(undef, M)
    ΥsinY = Array{Float64}(undef, M)
    ΥcosZ = Array{Float64}(undef, M)
    ΥsinZ = Array{Float64}(undef, M)
    
    #Initial measurement
    ρˣˣₖ₂[1] = gaugeStiffnessXX(k₂, ψ)
    ρˣˣₖ₃[1] = gaugeStiffnessXX(k₃, ψ)
    ρʸʸₖ₃[1] = gaugeStiffnessYY(k₃, ψ)
    ρʸʸₖ₁[1] = gaugeStiffnessYY(k₁, ψ)
    ρᶻᶻₖ₁[1] = gaugeStiffnessZZ(k₁, ψ)
    ρᶻᶻₖ₂[1] = gaugeStiffnessZZ(k₂, ψ)
    E_array[1] = E(ψ)
    ΥcosX[1], ΥcosY[1], ΥcosZ[1], ΥsinX[1], ΥsinY[1], ΥsinZ[1] = helicityModulusContribution(ψ)

    A_abs = 0.0 
    u⁺ = 0.0
    u⁻ = 0.0
    for x=1:L, y=1:L, z=1:L
        u⁺ += ψ.lattice[x,y,z].u⁺^2
        u⁻ += ψ.lattice[x,y,z].u⁻^2
        A_abs += ((ψ.lattice[x,y,z].A[1])^2 + (ψ.lattice[x,y,z].A[2])^2
                    + (ψ.lattice[x,y,z].A[3])^2 )
    end
    u⁺_array[1] = u⁺/(L^3 )
    u⁻_array[1] = u⁻/(L^3 )
    A_array[1] = A_abs/(L^3 )

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
        ΥcosX[m], ΥcosY[m], ΥcosZ[m], ΥsinX[m], ΥsinY[m], ΥsinZ[m] = helicityModulusContribution(ψ)

        A_abs = 0.0 
        u⁺ = 0.0
        u⁻ = 0.0
        for x=1:L, y=1:L, z=1:L
            u⁺ += ψ.lattice[x,y,z].u⁺
            u⁻ += ψ.lattice[x,y,z].u⁻
            A_abs += ((ψ.lattice[x,y,z].A[1])^2 + (ψ.lattice[x,y,z].A[2])^2
                    + (ψ.lattice[x,y,z].A[3])^2 )
        end
        u⁺_array[m] = u⁺/(L^3 )
        u⁻_array[m] = u⁻/(L^3 )
        A_array[m] = A_abs/(L^3 )
    end
    return (ρˣˣₖ₂, ρˣˣₖ₃, ρʸʸₖ₃, ρʸʸₖ₁, ρᶻᶻₖ₁, ρᶻᶻₖ₂, u⁺_array, u⁻_array, A_array, E_array,
            ΥcosX, ΥcosY, ΥcosZ, ΥsinX, ΥsinY, ΥsinZ)
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
        ΥcosX, ΥcosY, ΥcosZ, ΥsinX, ΥsinY, ΥsinZ) = gaugeStiffnessMeasure!(ψ_list[np], sim, M_min, Δt) 
   
    
    println("Measurements completed, collecting results from processes")
    #Gather results from workers
    for i=1:np-1
        (newρˣˣₖ₂, newρˣˣₖ₃, newρʸʸₖ₃, newρʸʸₖ₁, newρᶻᶻₖ₁, newρᶻᶻₖ₂, newu⁺, newu⁻, newA, newE,
            newΥcosX, newΥcosY, newΥcosZ, newΥsinX, newΥsinY, newΥsinZ) = fetch(futures[i])
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
        ΥcosX = vcat(ΥcosX, newΥcosX)
        ΥcosY = vcat(ΥcosY, newΥcosY)
        ΥcosZ = vcat(ΥcosZ, newΥcosZ)
        ΥsinX = vcat(ΥsinX, newΥsinX)
        ΥsinY = vcat(ΥsinY, newΥsinY)
        ΥsinZ = vcat(ΥsinZ, newΥsinZ)

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
    Υx = (mean(ΥcosX) - syst.β*mean(ΥsinX))*2/(syst.L)^2
    Υy = (mean(ΥcosY) - syst.β*mean(ΥsinY))*2/(syst.L)^2
    Υz = (mean(ΥcosZ) - syst.β*mean(ΥsinZ))*2/(syst.L₃)^2
    
    #Calculate the standard error 
    #TODO: Here we might also have to account for the autocorr time
    SEρˣˣₖ₂ = std(ρˣˣₖ₂)/sqrt(M)
    SEρˣˣₖ₃ = std(ρˣˣₖ₃)/sqrt(M)
    SEρʸʸₖ₃ = std(ρʸʸₖ₃)/sqrt(M)
    SEρʸʸₖ₁ = std(ρʸʸₖ₁)/sqrt(M)
    SEρᶻᶻₖ₁ = std(ρᶻᶻₖ₁)/sqrt(M)
    SEρᶻᶻₖ₂ = std(ρᶻᶻₖ₂)/sqrt(M)
    
    return (AVρˣˣₖ₂, SEρˣˣₖ₂, AVρˣˣₖ₃, SEρˣˣₖ₃, AVρʸʸₖ₃, SEρʸʸₖ₃, AVρʸʸₖ₁, SEρʸʸₖ₁,
            AVρᶻᶻₖ₁, SEρᶻᶻₖ₁, AVρᶻᶻₖ₂, SEρᶻᶻₖ₂, AVu⁺, AVu⁻, AVA, u⁺, AVE, Υx, Υy, Υz)
end

#--------------------------------------------------------------------------------------------------   
#Functionality to measure the helicity modulus 
function helicityModulusContribution(ψ::State)
    sumCosX = 0.0
    sumCosY = 0.0
    sumCosZ = 0.0
    sumSinX = 0.0
    sumSinY = 0.0
    sumSinZ = 0.0

    for i = 1:length(ψ.lattice)
        ϕ = ψ.lattice[i]
        ϕᵣ₊₁ = ψ.nb[i].ϕᵣ₊₁
        ϕᵣ₊₂ = ψ.nb[i].ϕᵣ₊₂
        ϕᵣ₊₃ = ψ.nb[i].ϕᵣ₊₃
        sumCosX += ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺)
        sumCosY += ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺)
        sumCosZ += ϕ.u⁺*ϕᵣ₊₃.u⁺*cos(ϕᵣ₊₃.θ⁺ - ϕ.θ⁺)
        sumSinX += ϕ.u⁺*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺)
        sumSinY += ϕ.u⁺*ϕᵣ₊₂.u⁺*sin(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺)
        sumSinZ += ϕ.u⁺*ϕᵣ₊₃.u⁺*sin(ϕᵣ₊₃.θ⁺ - ϕ.θ⁺)

    end
    return (sumCosX, sumCosY, sumCosZ, sumSinX^2, sumSinY^2, sumSinZ^2)
end

#--------------------------------------------------------------------------------------------------   
# Uses the helicityModulusContributions to calculate the helicity modulus vector
function helicityModulus(ψ::State)
    sum_cos_x, sum_cos_y, sum_cos_z, sum_sin_x2, sum_sin_y2, sum_sin_z2 = helicityModulusContribution(ψ)
    N = length(ψ.lattice)
    syst = ψ.consts
    Υ_x = 2/N*(sum_cos_x - 2*syst.β*sum_sin_x2);
    Υ_y = 2/N*(sum_cos_y - 2*syst.β*sum_sin_y2)
    Υ_z = 2/N*(sum_cos_z - 2*syst.β*sum_sin_z2)
    return [Υ_x, Υ_y, Υ_z]
end

#--------------------------------------------------------------------------------------------------   
# Assuming that the state consists of independent layers of 2D systems, calculate the helicity
# modulus Υ = (Υ_x + Υ_y)/2 for each layer, and return an average.
function helicityModulus2D(ψ::State)
    syst = ψ.consts; L = syst.L; L₃ = syst.L₃

    Υ = 0.0
    for z = 1:L₃

        sum_cos_x = 0.0; sum_cos_y = 0.0; sum_sin_x2 = 0.0; sum_sin_y2 = 0.0
        for v = 1:L, h = 1:L
            pos = [v, h, z]
            ϕ = ψ.lattice[pos...]
            ϕᵣ₊₁ = ψ.nb[pos...].ϕᵣ₊₁
            ϕᵣ₊₂ = ψ.nb[pos...].ϕᵣ₊₂
            
            sum_cos_x += ϕ.u⁺*ϕᵣ₊₁.u⁺*cos(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺)
            sum_cos_y += ϕ.u⁺*ϕᵣ₊₂.u⁺*cos(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺)
            sum_sin_x2 += ϕ.u⁺*ϕᵣ₊₁.u⁺*sin(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺)
            sum_sin_y2 += ϕ.u⁺*ϕᵣ₊₂.u⁺*sin(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺)
        end
    
        N = L^2
        Υ_x = 2/N*(sum_cos_x - 2*syst.β*(sum_sin_x2)^2)
        Υ_y = 2/N*(sum_cos_y - 2*syst.β*(sum_sin_y2)^2)

        Υ += (Υ_x + Υ_y)/2
    end

    return Υ/L₃
end

# --------------------------------------------------------------------------------------------------
# Calculates a list of helicity moduli Υ = (Υ_x + Υ_y)/2, one for each state, in parallel, and
# returns this list. If twoD == true, we assume the states consist of independent 2D systems.
function measureHelicityModulus(ψ_list::Array{State,1}; twoD=false)
    # Problem with memory management by usage of SharedArray and @sync @parallel for loop
    M = length(ψ_list)
#    Υ_list = SharedArray{Float64}(M)
    Υ_list = Array{Float64}(undef, M)
    if twoD
#        @sync @parallel for i = 1:M
        for i = 1:M
            Υ_list[i] = helicityModulus2D(ψ_list[i]::State)
        end
    else
#        @sync @parallel for i = 1:M
        for i = 1:M
            Υ_v = helicityModulus(ψ_list[i])
            Υ_list[i] = (Υ_v[1] + Υ_v[2])/2
        end
    end
#    return sdata(Υ_list)
    return Υ_list
end

###################################################################################################
#   Helicity modulus of the frozen 2-component superconductor, with only kinetic + maxwell terms
#
###################################################################################################

#--------------------------------------------------------------------------------------------------
# Takes a list of states and calculates the helicity modulus. Output returned as a three elements,
# one entry for each direction. This assumes that the hamiltonian is the frozen 2-component
# superconductor.
function helModFrozen2comp(ψ_list::Array{State,1})
    #Calculate the two average quantities needed to calculate the helicity modulus
    M = length(ψ_list)
    Υx_list = Array{Float64}(M)
    Υy_list = Array{Float64}(M)
    Υz_list = Array{Float64}(M)
    for i=1:M
        Υx_list[i], Υy_list[i], Υz_list[i] = helModF2compContributions(ψ_list[i])
    end
    return (Υx_list, Υy_list, Υz_list)
end

#--------------------------------------------------------------------------------------------------
# For a list of states, calculates the two average quantities required to calculate the helicity
# modulus, for each direction in the lattice.
function helModF2compContributions(ψ::State)
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    β = ψ.consts.β
    if L != L₃
        println("Warning from function helModF2compContributions: L and L₃ not equal")
    end
    
    sumcos_x = 0.0; sumcos_y = 0.0; sumcos_z = 0.0
    sumsin_x = 0.0; sumsin_y = 0.0; sumsin_z = 0.0         
    #Sum over lattice, calculating the contributions
    for h=1:L, v=1:L, z=1:L₃
        pos = [h,v,z]
        ϕ = ψ.lattice[pos...]
        ϕᵣ₊₁ = ψ.nb[pos...].ϕᵣ₊₁
        ϕᵣ₊₂ = ψ.nb[pos...].ϕᵣ₊₂
        ϕᵣ₊₃ = ψ.nb[pos...].ϕᵣ₊₃

        sumcos_x += 4*cos( (ϕᵣ₊₁.θ⁺ + ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - ϕ.θ⁻)/2.0 - ϕ.A[1] )*cos( (ϕᵣ₊₁.θ⁺ - ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ + ϕ.θ⁻)/2.0 )
        sumsin_x += 4*cos( (ϕᵣ₊₁.θ⁺ + ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ - ϕ.θ⁻)/2.0 - ϕ.A[1] )*sin( (ϕᵣ₊₁.θ⁺ - ϕᵣ₊₁.θ⁻ - ϕ.θ⁺ + ϕ.θ⁻)/2.0 )

        sumcos_y += 4*cos( (ϕᵣ₊₂.θ⁺ + ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - ϕ.θ⁻)/2.0 - ϕ.A[2] )*cos( (ϕᵣ₊₂.θ⁺ - ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ + ϕ.θ⁻)/2.0 )
        sumsin_y += 4*cos( (ϕᵣ₊₂.θ⁺ + ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ - ϕ.θ⁻)/2.0 - ϕ.A[2] )*sin( (ϕᵣ₊₂.θ⁺ - ϕᵣ₊₂.θ⁻ - ϕ.θ⁺ + ϕ.θ⁻)/2.0 )

        sumcos_z += 4*cos( (ϕᵣ₊₃.θ⁺ + ϕᵣ₊₃.θ⁻ - ϕ.θ⁺ - ϕ.θ⁻)/2.0 - ϕ.A[3] )*cos( (ϕᵣ₊₃.θ⁺ - ϕᵣ₊₃.θ⁻ - ϕ.θ⁺ + ϕ.θ⁻)/2.0 )
        sumsin_z += 4*cos( (ϕᵣ₊₃.θ⁺ + ϕᵣ₊₃.θ⁻ - ϕ.θ⁺ - ϕ.θ⁻)/2.0 - ϕ.A[3] )*sin( (ϕᵣ₊₃.θ⁺ - ϕᵣ₊₃.θ⁻ - ϕ.θ⁺ + ϕ.θ⁻)/2.0 )
        
    end
    norm = L^3
    Υx = (sumcos_x - β*sumsin_x^2)/norm
    Υy = (sumcos_y - β*sumsin_y^2)/norm
    Υz = (sumcos_z - β*sumsin_z^2)/norm
    return Υx, Υy, Υz
end

###################################################################################################
#       Order parameter(s) to measure the Chiral symmetry breaking
#
###################################################################################################


#Measuring when the phase difference is affected by Josephson ⟨cos(2(θ⁺-θ⁻))⟩
function josephsonResponse(ψ::State)
    L = ψ.consts.L
    sum = 0.0
    for h_pos=1:L, v_pos=1:L, z_pos=1:L
         ϕ = ψ.lattice[h_pos, v_pos, z_pos]
         sum += cos(2*(ϕ.θ⁺-ϕ.θ⁻))
    end
    j_response = sum/ψ.consts.L^3
    return j_response
end

function josephsonResponse(ψ_list::Array{State, 1})
    M = length(ψ_list)
    J_res = Array{Float64}(undef, M)
    for m=1:M
        J_res = josephsonResponse(ψ_list[m])
    end
    return J_res
end
    

function measurementSeriesE!(ψ::State, M::Int64, Δt::Int64, sim::Controls)
    L = ψ.consts.L
    E_array=Array{Float64, 1}(undef, M)
    θ_array = Array{Float64, 1}(undef, M)
    ρˣˣₖ₂ = Array{Float64, 1}(undef, M)
    ρˣˣₖ₃ = Array{Float64, 1}(undef, M)
    ρʸʸₖ₁ = Array{Float64, 1}(undef, M)
    ρʸʸₖ₃ = Array{Float64, 1}(undef, M)
    ρᶻᶻₖ₁ = Array{Float64, 1}(undef, M)
    ρᶻᶻₖ₂ = Array{Float64, 1}(undef, M)
    k₁ = [2π/L, 0.0, 0.0]
    k₂ = [0.0, 2π/L, 0.0]
    k₃ = [0.0, 0.0, 2π/L]
    
    if ψ.PBC
        E_array[1] = E(ψ)
        θ_array[1] = josephsonResponse(ψ)
         ρˣˣₖ₂[1], _, ρᶻᶻₖ₂[1] = gaugeStiffness(k₂, ψ)
        ρˣˣₖ₃[1], ρʸʸₖ₃[1], _ = gaugeStiffness(k₃, ψ)
        _, ρʸʸₖ₁[1], ρᶻᶻₖ₁[1] = gaugeStiffness(k₁, ψ)
    else
        E_array[1] = E_APBC(ψ)
        θ_array[1] = josephsonResponse(ψ)
        ρˣˣₖ₂[1], _, ρᶻᶻₖ₂[1] = gaugeStiffness(k₂, ψ)
        ρˣˣₖ₃[1], ρʸʸₖ₃[1], _ = gaugeStiffness(k₃, ψ)
        _, ρʸʸₖ₁[1], ρᶻᶻₖ₁[1] = gaugeStiffness(k₁, ψ)
    end
                
    for i=2:M
        for j=1:Δt
            mcSweep!(ψ, sim)
        end
        if ψ.PBC
            E_array[i] = E(ψ)
            θ_array[i] = josephsonResponse(ψ)
            ρˣˣₖ₂[i], _, ρᶻᶻₖ₂[i] = gaugeStiffness(k₂, ψ)
            ρˣˣₖ₃[i], ρʸʸₖ₃[i], _ = gaugeStiffness(k₃, ψ)
            _, ρʸʸₖ₁[i], ρᶻᶻₖ₁[i] = gaugeStiffness(k₁, ψ)
        else    
            E_array[i] = E_APBC(ψ)
            θ_array[i] = josephsonResponse(ψ)
            ρˣˣₖ₂[i], _, ρᶻᶻₖ₂[i] = gaugeStiffness(k₂, ψ)
            ρˣˣₖ₃[i], ρʸʸₖ₃[i], _ = gaugeStiffness(k₃, ψ)
            _, ρʸʸₖ₁[i], ρᶻᶻₖ₁[i] = gaugeStiffness(k₁, ψ)
        end
    end
    return E_array, θ_array, ρˣˣₖ₂, ρˣˣₖ₃, ρʸʸₖ₁, ρʸʸₖ₃, ρᶻᶻₖ₁, ρᶻᶻₖ₂
end

function measurementSeriesEParalell!(ψ_list::Array{State,1}, M::Int64, Δt::Int64, sim::Controls)
    N_states = length(ψ_list)
    m = Int64(floor(M/N_states))
    E_matrix = Array{Float64, 2}(undef, N_states-1, m)
    θ_matrix = Array{Float64, 2}(undef, N_states-1, m)
    ρˣˣₖ₂ = Array{Float64, 2}(undef, N_states-1, m)
    ρˣˣₖ₃ = Array{Float64, 2}(undef, N_states-1, m)
    ρʸʸₖ₁ = Array{Float64, 2}(undef, N_states-1, m)
    ρʸʸₖ₃ = Array{Float64, 2}(undef, N_states-1, m)
    ρᶻᶻₖ₁ = Array{Float64, 2}(undef, N_states-1, m)
    ρᶻᶻₖ₂ = Array{Float64, 2}(undef, N_states-1, m)

    E_future_list = [Future() for i=1:N_states-1]
    for i=1:N_states-1
        E_future_list[i] = @spawn measurementSeriesE!(ψ_list[i], m, Δt, sim)
    end
    (E_last, θ_last, ρˣˣₖ₂_last, ρˣˣₖ₃_last, ρʸʸₖ₁_last, ρʸʸₖ₃_last, 
        ρᶻᶻₖ₁_last, ρᶻᶻₖ₂_last) = measurementSeriesE!(ψ_list[N_states], m, Δt, sim)
    for i=1:N_states-1
        (E_matrix[i,:], θ_matrix[i,:], ρˣˣₖ₂[i,:], ρˣˣₖ₃[i,:], ρʸʸₖ₁[i,:], ρʸʸₖ₃[i,:], 
            ρᶻᶻₖ₁[i,:], ρᶻᶻₖ₂[i,:]) = fetch(E_future_list[i])
    end
    for i=1:N_states-1
        E_last = vcat(E_last, E_matrix[i,:])
        θ_last = vcat(θ_last, θ_matrix[i,:])
        ρˣˣₖ₂_last = vcat(ρˣˣₖ₂_last, ρˣˣₖ₂[i,:])
        ρˣˣₖ₃_last = vcat(ρˣˣₖ₃_last, ρˣˣₖ₃[i,:])
        ρʸʸₖ₁_last = vcat(ρʸʸₖ₁_last, ρʸʸₖ₁[i,:])
        ρʸʸₖ₃_last = vcat(ρʸʸₖ₃_last, ρʸʸₖ₃[i,:]) 
        ρᶻᶻₖ₁_last = vcat(ρᶻᶻₖ₁_last, ρᶻᶻₖ₁[i,:]) 
        ρᶻᶻₖ₂_last = vcat(ρᶻᶻₖ₂_last, ρᶻᶻₖ₂[i,:])

    end
    return E_last, θ_last, ρˣˣₖ₂_last, ρˣˣₖ₃_last, ρʸʸₖ₁_last, ρʸʸₖ₃_last, ρᶻᶻₖ₁_last, ρᶻᶻₖ₂_last 
end
function lineTensionFromHighT(Eₛ::Array{Float64,1}, β_array::Array{Float64,1}, N::Int64, Δβ::Float64, L::Int64)
    lineTension = Array{Float64, 1}(undef, N)
    lineTension[1] = 0.0
    for i=2:N
        lineTension[i] = lineTension[i-1] + Δβ*(Eₛ[i]+Eₛ[i-1])/2.0
    end
    
    for i=1:N
        lineTension[i] = lineTension[i]/(β_array[i]*L^2)
    end
    return lineTension
end

function lineTensionFromHighT_Norm(Eₛ::Array{Float64,1}, β_array::Array{Float64,1}, N::Int64, Δβ::Float64, L::Int64)
    lineTension = Array{Float64, 1}(undef, N)
    lineTension[1] = 0.0
    for i=2:N
        lineTension[i] = lineTension[i-1] + Δβ*(Eₛ[i]+Eₛ[i-1])/2.0
    end
    
    for i=1:N
        lineTension[i] = lineTension[i]/L^2
    end
    return lineTension
end
