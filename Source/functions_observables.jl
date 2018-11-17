####################################################################################################################
#                            Planar Vorticity
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# New version of nᵣ based on suggestion by Troels. We are drawing the gauge-invariant phase
# difference back to [-π, π) instead of [0, 2π) and also adding 2πf so that we are measuring
# n instead of (n-f) which we would do by just doing the gauge-invariant phase difference.
function drawback{T<:Real}(x::T)
    return mod2pi(x+π)-π
end
function nᵣ(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
    vort_θ⁺ = (drawback(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ - ϕ.A[1]) + drawback(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₁.θ⁺ - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos))
        - drawback(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₂.θ⁺  - ϕᵣ₊₂.A[1])
        - drawback(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺  - (ϕ.A[2] + two_pi*c.f*(h_pos-1)))+two_pi*c.f)
    vort_θ⁻ = (drawback(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ - ϕ.A[1]) + drawback(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₁.θ⁻ - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos))
        - drawback(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₂.θ⁻  - ϕᵣ₊₂.A[1])
        - drawback(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻  - (ϕ.A[2] + two_pi*c.f*(h_pos-1)))+two_pi*c.f)
    return vort_θ⁺, vort_θ⁻
end
function nᵣNoA(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
    vort_θ⁺ = (mod(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺, two_pi) + mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₁.θ⁺, two_pi)
        - mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₂.θ⁺, two_pi)
        - mod(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺, two_pi))
    vort_θ⁻ = (mod(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻, two_pi) + mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₁.θ⁻, two_pi)
        - mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₂.θ⁻, two_pi) 
        - mod(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻, two_pi))
    return vort_θ⁺, vort_θ⁻
end
function symmetrizedVorticity(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, ϕᵣ₊₂₂::LatticeSite, ϕᵣ₋₁₊₂::LatticeSite, h_pos::Int64)
    Aᵣ = c.f*two_pi*(h_pos-1)
    vort_θ⁺ = (drawback(ϕᵣ₊₁₊₂.θ⁺ - ϕ.θ⁺ - (ϕ.A[2]+Aᵣ+ϕᵣ₊₂.A[1]) ) + drawback(ϕᵣ₊₂₂.θ⁺ - ϕᵣ₊₁₊₂.θ⁺ - (-ϕᵣ₊₂.A[1]+ϕᵣ₊₂.A[2]+Aᵣ))
               - drawback(ϕᵣ₊₂₂.θ⁺ - ϕᵣ₋₁₊₂.θ⁺ - (ϕᵣ₋₁₊₂.A[1]+ϕᵣ₊₂.A[2]+Aᵣ)) - drawback(ϕᵣ₋₁₊₂.θ⁺ - ϕ.θ⁺ - (ϕ.A[2]+Aᵣ-ϕᵣ₋₁₊₂.A[1])))
    vort_θ⁻ = (drawback(ϕᵣ₊₁₊₂.θ⁻ - ϕ.θ⁻ - (ϕ.A[2]+Aᵣ+ϕᵣ₊₂.A[1]) ) + drawback(ϕᵣ₊₂₂.θ⁻ - ϕᵣ₊₁₊₂.θ⁻ - (-ϕᵣ₊₂.A[1]+ϕᵣ₊₂.A[2]+Aᵣ))
               - drawback(ϕᵣ₊₂₂.θ⁻ - ϕᵣ₋₁₊₂.θ⁻ - (ϕᵣ₋₁₊₂.A[1]+ϕᵣ₊₂.A[2]+Aᵣ)) - drawback(ϕᵣ₋₁₊₂.θ⁻ - ϕ.θ⁻ - (ϕ.A[2]+Aᵣ-ϕᵣ₋₁₊₂.A[1])))
    return vort_θ⁺, vort_θ⁻
end


# -----------------------------------------------------------------------------------------------------------
# Assuming we have a state ψ, we want to find the lattice of vortices.
function vortexSnapshot(ψ::State)
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    V⁺ = zeros(L,L,L₃)
    V⁻ = zeros(L,L,L₃)
    
    # Sum over the lattice
    for z_pos = 1:L₃, h_pos = 1:L, v_pos = 1:L
        ϕ = ψ.lattice[v_pos,h_pos,z_pos]
#        ϕᵣ₊₁₊₂ = ψ.nnb[v_pos,h_pos,z_pos].ϕᵣ₊₁₊₂
#        ϕᵣ₋₁₊₂ = ψ.nnb[v_pos,h_pos,z_pos].ϕᵣ₋₁₊₂
#        ϕᵣ₊₂₂ = ψ.nnnb[v_pos,h_pos,z_pos].ϕᵣ₊₂₂
        ϕᵣ₊₁ = ψ.nb[v_pos,h_pos,z_pos].ϕᵣ₊₁
        ϕᵣ₊₂ = ψ.nb[v_pos,h_pos,z_pos].ϕᵣ₊₂
        ϕᵣ₊₁₊₂ = ψ.nnb[v_pos,h_pos,z_pos].ϕᵣ₊₁₊₂
        (V⁺[v_pos,h_pos,z_pos], V⁻[v_pos, h_pos, z_pos]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)
#        V⁺[v_pos,h_pos,z_pos], V⁻[v_pos,h_pos,z_pos] = symmetrizedVorticity(ψ.consts, ϕ, ϕᵣ₊₁₊₂, ϕᵣ₊₂₂, ϕᵣ₋₁₊₂)
    end
    
    return (V⁺, V⁻)
end

# -----------------------------------------------------------------------------------------------------------
# We combine the two component vortex lattices into one lattice where
# 0: (-1, -1), 1: (-1, 0), 2: (-1, 1)
# 3: (0, -1), 4: (0, 0), 5: (0, 1)
# 6: (1, -1), 7: (1, 0), 8: (1, 1)
function combineVortexLattices{T<:Real}(vortex_matrix⁺::Array{T, 2}, vortex_matrix⁻::Array{T,2})
    L = size(vortex_matrix⁺,1)
    A = [-1 for x=1:L, y=1:L]
    for v_pos = 1:L, h_pos = 1:L
        if isapprox(vortex_matrix⁺[v_pos, h_pos], -1.0, atol=0.3)
            if isapprox(vortex_matrix⁻[v_pos, h_pos], -1.0, atol=0.3)
                A[v_pos, h_pos] = 0
            elseif isapprox(vortex_matrix⁻[v_pos, h_pos], 0.0, atol=0.3)
                A[v_pos, h_pos] = 1
            elseif isapprox(vortex_matrix⁻[v_pos, h_pos], 1.0, atol=0.3)
                A[v_pos, h_pos] = 2
            end
        elseif isapprox(vortex_matrix⁺[v_pos, h_pos], 0.0, atol=0.3)
            if isapprox(vortex_matrix⁻[v_pos, h_pos], -1.0, atol=0.3)
                A[v_pos, h_pos] = 3
            elseif isapprox(vortex_matrix⁻[v_pos, h_pos], 0.0, atol=0.3)
                A[v_pos, h_pos] = 4
            elseif isapprox(vortex_matrix⁻[v_pos, h_pos], 1.0, atol=0.3)
                A[v_pos, h_pos] = 5
            end
        elseif isapprox(vortex_matrix⁺[v_pos, h_pos], 1.0, atol=0.3)
            if isapprox(vortex_matrix⁻[v_pos, h_pos], -1.0, atol=0.3)
                A[v_pos, h_pos] = 6
            elseif isapprox(vortex_matrix⁻[v_pos, h_pos], 0.0, atol=0.3)
                A[v_pos, h_pos] = 7
            elseif isapprox(vortex_matrix⁻[v_pos, h_pos], 1.0, atol=0.3)
                A[v_pos, h_pos] = 8
            end
        end
    end
    return A
end

# --------------------------------------------------------------------------------------------------
# Projects 3D vorticity space down on 2D through an averaging over the z-dimension.
function avgVort{T<:Real}(V::Array{T,3})
    L = size(V,1)
    L₃ = size(V,3)
    avg_V = zeros(L,L)
    for h_pos = 1:L, v_pos = 1:L
        for z_pos = 1:L₃
            avg_V[v_pos,h_pos] += V[v_pos,h_pos,z_pos]
        end
        avg_V[v_pos,h_pos] /= L₃
    end
    return avg_V
end

####################################################################################################################
#                            Structure function
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# V⁺ and V⁻ are LxL matrices of plaquettes containing the vorticities nᵣ of the two components averaged over
# the z-direction.
# Returns the structure function at k of both vorticities.
function structureFunction{T<:Real, I<:Real}(k::Array{T,1}, ψ::State, V⁺::Array{I,2}, V⁻::Array{I,2})
    sum⁺ = Complex(0)
    sum⁻ = Complex(0)
    L = ψ.consts.L
    origo_l = Int(floor((L-1)/2))
    origo_r = L-1-origo_l # For a lattice where L is odd, origo will be in the center.
    
    # Sum over all the plaquettes in the whole lattice
    for h_pos = 1:L
        for v_pos = 1:L
            r = [h_pos-1-origo_l, origo_r+1-v_pos]
            sum⁺ += V⁺[v_pos,h_pos]*exp(im*(k⋅r))
            sum⁻ += V⁻[v_pos,h_pos]*exp(im*(k⋅r))
        end
    end
    
    return (abs2(sum⁺),abs2(sum⁻))
end
# If V⁺ and V⁻ are not known, they have to be calculated.
function structureFunction{T<:Real}(k::Array{T,1}, ψ::State)
    sum⁺ = Complex(0)
    sum⁻ = Complex(0)
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    
    # Sum over lattice
    for h_pos = 1:L, v_pos = 1:L
        r = [h_pos-1, L-v_pos]		# For r we assume origo is in position [L,1] of the lattice. 
              # Note that r is the same as pos (found previously) with y-axis flipped and -1 in each direction.
              # Additionally we define it such that we get the usual r = [x,y] order of dimensions.

        vort_θ⁺ = 0.0; vort_θ⁻ = 0.0
        for z_pos = 1:L₃
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.nb[v_pos,h_pos].ϕᵣ₊₁
            ϕᵣ₊₂ = ψ.nb[v_pos,h_pos].ϕᵣ₊₂
            ϕᵣ₊₁₊₂ = ψ.nnb[v_pos,h_pos].ϕᵣ₊₁₊₂
            n⁺, n⁻ = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)
            vort_θ⁺ += n⁺
            vort_θ⁻ += n⁻
        end
        vort_θ⁺ /= L₃
        vort_θ⁻ /= L₃
        sum⁺ += vort_θ⁺*exp(im*(k⋅r))
        sum⁻ += vort_θ⁻*exp(im*(k⋅r))
    end
    
    return (abs2(sum⁺),abs2(sum⁻))
end

####################################################################################################################
#                            Thermal averages
#
####################################################################################################################

# --------------------------------------------------------------------------------------------------
# Take in a matrix of k-values and calculate both the vorticity of θ⁺ and θ⁻, as well as an average over the
# real space vortex lattice.
# Assumes input state is at equilibrium.
function structureFunctionVortexLatticeAvg!{T<:Real}(ks::Array{Array{T, 1}, 2}, 
        ψ::State, sim::Controls, M::Int64, Δt::Int64)
    syst = ψ.consts
    L = syst.L
    
    # Setup structure factor storage
    Lky = size(ks, 1)
    Lkx = size(ks, 2)
    S⁺ = [zeros(Lky,Lkx) for i=1:M] # Series of matrices, one matrix for each measurement.
    S⁻ = [zeros(Lky,Lkx) for i=1:M]
    Sm⁺ = zeros(Lky,Lkx)
    Sm⁻ = zeros(Lky,Lkx)
    s_norm_inv = 1/(L^2*syst.f*two_pi)^2
    
    # Setup vortex lattice storage
    V⁺ = [zeros(L,L) for i=1:M]    # Matrix containing the series of measurements for each position
    V⁻ = [zeros(L,L) for i=1:M]
    VSm⁺ = zeros(L,L)
    VSm⁻ = zeros(L,L)
    avV⁺ = zeros(L,L)
    avV⁻ = zeros(L,L)
    
    println("Making $(M) measurements using $(M*Δt) MCS over a $(Lkx)×$(Lky) matrix of ks.")
    # Loop over M measurements
    for m = 1:M
        print("Measurement progress: $(Int(round(m/M*100,0)))% \r")
        flush(STDOUT)
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ, sim)
        end
        
        # Find n_z(r) of the lattice.
        (V⁺[m], V⁻[m]) = vortexSnapshot(ψ)
        
        # Calculate average of vorticity second moment. 
        for y = 1:L, x = 1:L
            VSm⁺[y,x] += V⁺[m][y,x]^2
            VSm⁻[y,x] += V⁻[m][y,x]^2
        end
        
        # Find structure factor. 
        for y = 1:Lky, x = 1:Lkx
            (S⁺[m][y,x], S⁻[m][y,x]) = s_norm_inv.*structureFunction(ks[y,x], ψ, V⁺[m], V⁻[m])
            Sm⁺[y,x] += S⁺[m][y,x]^2
            Sm⁻[y,x] += S⁻[m][y,x]^2
        end
        
        V⁺[m] = V⁺[m]./two_pi
        V⁻[m] = V⁻[m]./two_pi
    end
    
    # Error calculation of average vorticity
    τ_V⁺ = [autocorrTime([V⁺[m][y,x] for m=1:M], 5.0) for y=1:L, x=1:L]
    τ_V⁻ = [autocorrTime([V⁻[m][y,x] for m=1:M], 5.0) for y=1:L, x=1:L]
    errV⁺ = zeros(L,L)
    errV⁻ = zeros(L,L)
    
    avV⁺ = mean(V⁺)
    avV⁻ = mean(V⁻)
    
    for y=1:L, x=1:L
        VSm⁺[y,x] /= M
        VSm⁻[y,x] /= M
        errV⁺[y,x] = (1+2*τ_V⁺[y,x])*(VSm⁺[y,x] - avV⁺[y,x]^2)/(M-1)
        errV⁻[y,x] = (1+2*τ_V⁻[y,x])*(VSm⁻[y,x] - avV⁻[y,x]^2)/(M-1)
    end
    
    # Error calculation of structure factor.
    avS⁺ = mean(S⁺)
    avS⁻ = mean(S⁻)
    τ_S⁺ = [autocorrTime([S⁺[m][y,x] for m=1:M], 5.0) for y=1:Lky, x=1:Lkx]
    τ_S⁻ = [autocorrTime([S⁻[m][y,x] for m=1:M], 5.0) for y=1:Lky, x=1:Lkx]
    errS⁺ = [0.0 for y=1:Lky, x=1:Lkx]
    errS⁻ = [0.0 for y=1:Lky, x=1:Lkx]
    
    for y=1:Lky, x=1:Lkx
        Sm⁺[y,x] /= M
        Sm⁻[y,x] /= M
        errS⁺[y,x] = (1+2*τ_S⁺[y,x])*(Sm⁺[y,x] - avS⁺[y,x]^2)/(M-1)
        errS⁻[y,x] = (1+2*τ_S⁻[y,x])*(Sm⁻[y,x] - avS⁻[y,x]^2)/(M-1)
    end
    
    # Finding max values over the matrices.
    max_S⁺= maximum(avS⁺)
    max_S⁻ = maximum(avS⁻)
    max_err_S⁺ = maximum(errS⁺)
    max_err_S⁻ = maximum(errS⁻)
    max_τ_S⁺ = maximum(τ_S⁺)
    max_τ_S⁻ = maximum(τ_S⁻)
    println("\nMax (S⁺, S⁻)\n($(max_S⁺), $(max_S⁻))")
    println("Max δ(S⁺, S⁻)\n($(max_err_S⁺), $(max_err_S⁻))")
    println("Max correlation time\n($(max_τ_S⁺), $(max_τ_S⁻))")
    
    return (avV⁺, errV⁺, V⁺, avV⁻, errV⁻, V⁻, avS⁺, errS⁺, S⁺, avS⁻, errS⁻, S⁻)
end


# --------------------------------------------------------------------------------------------------
# Makes M measurements of both the vortex lattice and the structure factor.
function sfvlaMeasure!{T<:Real}(ks::Array{Array{T, 1}, 2}, ψ::State, sim::Controls, M::Int64, Δt::Int64)
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    L_k = size(ks, 1)
    # Setup measurement storage
    S⁺ = [zeros(L_k,L_k) for i=1:M]
    S⁻ = [zeros(L_k,L_k) for i=1:M]
    V⁺ = [zeros(L,L,L₃) for i=1:M]
    V⁻ = [zeros(L,L,L₃) for i=1:M]
    proj_V⁺ = [zeros(L,L) for i=1:M]
    proj_V⁻ = [zeros(L,L) for i=1:M]
    
    # The first measurement is of the initial ψ
    (V⁺[1], V⁻[1]) = vortexSnapshot(ψ)
    proj_V⁺[1] = avgVort(V⁺[1]); proj_V⁻[1] = avgVort(V⁻[1])
    
    # Then we use this to measure the structure factor
    for x=1:L_k, y=1:L_k
        (S⁺[1][y,x], S⁻[1][y,x]) = structureFunction(ks[y,x], ψ, proj_V⁺[1], proj_V⁻[1])
    end
    
    # For each measurement
    for m = 2:M
#        print("Measurement progress: $(Int(round(m/M*100,0)))% \r")
#        flush(STDOUT)
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ, sim)
        end
        
        # Find n_z(r) of the lattice.
        (V⁺[m], V⁻[m]) = vortexSnapshot(ψ)
        # Project down on a 2D plane
        proj_V⁺[m] = avgVort(V⁺[m]); proj_V⁻[m] = avgVort(V⁻[m])
        # Find structure factor. 
        for x=1:L_k, y=1:L_k
            (S⁺[m][y,x], S⁻[m][y,x]) = structureFunction(ks[y,x], ψ, proj_V⁺[m], proj_V⁻[m])
        end
    end
    
    # After the loop, we should have filled up the M measurements for the matrices.
    return (S⁺, S⁻, proj_V⁺, proj_V⁻)
end

# Same as above, but now prints progress to STDOUT
function sfvlaMeasure!{T<:Real}(ks::Array{Array{T, 1}, 2}, ψ::State, sim::Controls, M::Int64,
        Δt::Int64, option::AbstractString)
	PROG_NUM = 10		# Adjusts the number of times progress is reported while measuring.
    L = ψ.consts.L
    L₃ = ψ.consts.L₃
    L_k = size(ks, 1)
    # Setup measurement storage
    S⁺ = [zeros(L_k,L_k) for i=1:M]
    S⁻ = [zeros(L_k,L_k) for i=1:M]
    V⁺ = [zeros(L,L,L₃) for i=1:M]
    V⁻ = [zeros(L,L,L₃) for i=1:M]
    proj_V⁺ = [zeros(L,L) for i=1:M]
    proj_V⁻ = [zeros(L,L) for i=1:M]
    
    # The first measurement is of the initial ψ
    (V⁺[1], V⁻[1]) = vortexSnapshot(ψ)
    proj_V⁺[1] = avgVort(V⁺[1]); proj_V⁻[1] = avgVort(V⁻[1])
    
    # Then we use this to measure the structure factor
    for x=1:L_k, y=1:L_k
        (S⁺[1][y,x], S⁻[1][y,x]) = structureFunction(ks[y,x], ψ, proj_V⁺[1], proj_V⁻[1])
    end
    
    # For each measurement
	prog_int = floor(Int64, M/PROG_NUM)
    for m = 2:M
		if m % prog_int == 0
			println("Measurement progress: $(Int(round(m/M*100,0)))%")
        	flush(STDOUT)
		end
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ, sim)
        end
        
        # Find n_z(r) of the lattice.
        (V⁺[m], V⁻[m]) = vortexSnapshot(ψ)
        proj_V⁺[m] = avgVort(V⁺[m]); proj_V⁻[m] = avgVort(V⁻[m])
        
        # Find structure factor. 
        for x=1:L_k, y=1:L_k
            (S⁺[m][y,x], S⁻[m][y,x]) = structureFunction(ks[y,x], ψ, proj_V⁺[m], proj_V⁻[m])
        end
    end
    
    # After the loop, we should have filled up the M measurements for the matrices.
    return (S⁺, S⁻, proj_V⁺, proj_V⁻)
end

# --------------------------------------------------------------------------------------------------
# Given np+1 uncorrelated states in ψ_list we use these to make M measurements of the vortex lattice by splitting
# the M measurements on the np workers as well as the master process. In this version we continuously discard
# the measurements and only save the averages and second moments.
function parallelSFVLA!{T<:Real}(ks::Array{Array{T, 1}, 2}, 
        ψ_list::Array{State,1}, sim::Controls, M::Int64, Δt::Int64)
    syst = ψ_list[1].consts
    L = syst.L
    
    # Checking that the k matrix has equal dimensions
    L_k = size(ks, 1)
    size(ks, 2) == L_k || throw(DomainError())
    
    s_norm_inv = 1/(L^2*syst.f*two_pi)^2
    v_norm_inv = 1/two_pi
    
    # Splitting the problem into np sub-problems.
    np = nprocs()
    # Minimum amount of work pr. process
    M_min = Int(floor(M/np))
    # Number of workers doing +1 extra work
    nw = M%np
    
    # Make sure that we have enough states
    length(ψ_list) >= np || throw(error("ERROR: Not enough states in list"))
    
    # Setup worker futures
    futures = [Future() for i=1:(np-1)]
    
    println("Starting $(M) measurements on $(np) processes doing max $(M_min + Int(ceil(nw/np))) measurements each
on a $(L)×$(L) system, corresponding to $((M_min+ceil(Int64, nw/np))*Δt) MCS pr. process")
    
    # Start +1 workers
    for i = 1:nw
        futures[i] = @spawn sfvlaMeasure!(ks, ψ_list[i], sim, M_min+1, Δt)
    end
    # Start remaining workers
    for i = 1:np-nw-1
        futures[nw+i] = @spawn sfvlaMeasure!(ks, ψ_list[nw+i], sim, M_min, Δt)
    end
    # Make the master process work as well
    (S⁺, S⁻, V⁺, V⁻) = sfvlaMeasure!(ks, ψ_list[np], sim, M_min, Δt, "-v")
    
    println("Measurements done, collecting parallell results.")
    # Collect results
    for i = 1:np-1
        (new_S⁺, new_S⁻, new_V⁺, new_V⁻) = fetch(futures[i])
        S⁺ = vcat(S⁺, new_S⁺)
        S⁻ = vcat(S⁻, new_S⁻)
        V⁺ = vcat(V⁺, new_V⁺)
        V⁻ = vcat(V⁻, new_V⁻)
    end
    @test length(S⁺) == M
    
    println("Parallell measurements done. Processing.")
    
    # Normalizing results
    V⁺ = v_norm_inv.*V⁺
    V⁻ = v_norm_inv.*V⁻
    S⁺ = s_norm_inv.*S⁺
    S⁻ = s_norm_inv.*S⁻
    
    # Calculate averages and second moments
    av_S⁺ = mean(S⁺)
    av_S⁻ = mean(S⁻)
    av_V⁺ = mean(V⁺)
    av_V⁻ = mean(V⁻)
    sm_S⁺ = zeros(L_k, L_k)
    sm_S⁻ = zeros(L_k, L_k)
    sm_V⁺ = zeros(L,L)
    sm_V⁻ = zeros(L,L)
    for m = 1:M
        sm_S⁺ += S⁺[m].^2
        sm_S⁻ += S⁻[m].^2
        sm_V⁺ += V⁺[m].^2
        sm_V⁻ += V⁻[m].^2
    end
    sm_S⁺ = sm_S⁺./M
    sm_S⁻ = sm_S⁻./M;
    sm_V⁺ = sm_V⁺./M
    sm_V⁻ = sm_V⁻./M
    
    # Error calculation of average vorticity
    τ_V⁺ = [1/ess_factor_estimate([V⁺[m][y,x] for m=1:M])[1] for y=1:L, x=1:L]
    τ_V⁻ = [1/ess_factor_estimate([V⁻[m][y,x] for m=1:M])[1] for y=1:L, x=1:L]
    
    err_V⁺ = (1+2.*τ_V⁺).*(sm_V⁺ - av_V⁺.^2)./(M-1)
    err_V⁻ = (1+2.*τ_V⁻).*(sm_V⁻ - av_V⁻.^2)./(M-1)
    
    # Error calculation of structure factor.
    τ_S⁺ = [1/ess_factor_estimate([S⁺[m][y,x] for m=1:M])[1] for y=1:L_k, x=1:L_k]
    τ_S⁻ = [1/ess_factor_estimate([S⁻[m][y,x] for m=1:M])[1] for y=1:L_k, x=1:L_k]
    
    err_S⁺ = (1+2.*τ_S⁺).*(sm_S⁺ - av_S⁺.^2)./(M-1)
    err_S⁻ = (1+2.*τ_S⁻).*(sm_S⁻ - av_S⁻.^2)./(M-1)

	# Sum of all vorticities
	println("\nSum of vorticity of random snapshot:\nV⁺: \t$(sum(V⁺[rand(1:M)]))
V⁻: \t$(sum(V⁻[rand(1:M)]))")
    
    
    # Finding max values over the matrices.
    max_S⁺= maximum(av_S⁺)
    max_S⁻ = maximum(av_S⁻)
    max_err_S⁺ = maximum(err_S⁺)
    max_err_S⁻ = maximum(err_S⁻)
    max_τ_S⁺ = maximum(τ_S⁺)
    max_τ_S⁻ = maximum(τ_S⁻)
    println("\nMax (S⁺, S⁻)\n($(max_S⁺), $(max_S⁻))")
    println("Max δ(S⁺, S⁻)\n($(max_err_S⁺), $(max_err_S⁻))")
    println("Max correlation time\n($(max_τ_S⁺), $(max_τ_S⁻))")
    
    return (av_V⁺, err_V⁺, V⁺, av_V⁻, err_V⁻, V⁻, av_S⁺, err_S⁺, S⁺, av_S⁻, err_S⁻, S⁻)
end



####################################################################################################################
#                            Amplitude measurements
#
####################################################################################################################

function avgAmpMeasure!(ψ::State, sim::Controls, M::Int64, Δt::Int64; visible=false)
    
    PROG_NUM = 10
    u⁺_list = Array{Float64, 1}(M); u⁻_list = Array{Float64, 1}(M)
    u⁺_list[1], u⁻_list[1] = meanAmplitudes(ψ)
    
    prog_int = floor(Int64, M/PROG_NUM)
    if visible
        for m = 2:M
            nMCS(ψ, sim, Δt)
            u⁺_list[m], u⁻_list[m] = meanAmplitudes(ψ)
            if m % prog_int == 0
                println("Measurement progress: $(Int(round(m/M*100,0)))%")
                flush(STDOUT)
            end
        end
    else # If not visible
        for m = 2:M
            nMCS(ψ, sim, Δt)
            u⁺_list[m], u⁻_list[m] = meanAmplitudes(ψ)
        end
    end
    
    return u⁺_list, u⁻_list
end

# --------------------------------------------------------------------------------------------------
# Given np+1 uncorrelated states in ψ_list we use these to make M measurements of the vortex lattice by splitting
# the M measurements on the np workers as well as the master process. In this version we continuously discard
# the measurements and only save the averages and second moments.
function averageAmplitudes!(ψ_list::Array{State,1}, sim::Controls, M::Int64, Δt::Int64)
    syst = ψ_list[1].consts
    L = syst.L
    
    # Setup storage for measurements of average amplitudes
    u⁺_list = Array{Float64, 1}(M); u⁻_list = Array{Float64, 1}(M);
    
    # Splitting the problem into np sub-problems.
    np = nprocs()
    # Minimum amount of work pr. process
    M_min = Int(floor(M/np))
    # Number of workers doing +1 extra work
    nw = M%np
    
    # Make sure that we have enough states
    length(ψ_list) >= np || throw(error("ERROR: Not enough states in list"))
    
    # Setup worker futures
    futures = [Future() for i=1:(np-1)]
    
    println("Starting $(M) measurements on $(np) processes doing max $(M_min + Int(ceil(nw/np))) measurements each
on a $(L)×$(L) system, corresponding to $((M_min+ceil(Int64, nw/np))*Δt) MCS pr. process")
    
    # Start +1 workers
    for i = 1:nw
        futures[i] = @spawn avgAmpMeasure!(ψ_list[i], sim, M_min+1, Δt)
    end
    # Start remaining workers
    for i = 1:np-nw-1
        futures[nw+i] = @spawn avgAmpMeasure!(ψ_list[nw+i], sim, M_min, Δt)
    end
    # Make the master process work as well
    
    u⁺_list[1:M_min], u⁻_list[1:M_min] = avgAmpMeasure!(ψ_list[np], sim, M_min, Δt; visible=true)
    
    println("Measurements done, collecting parallell results.")
    # Collect results
    for i = 1:nw
        int = (i-1)*(M_min+1)+1+M_min:i*(M_min+1)+M_min
        u⁺_list[int], u⁻_list[int] = fetch(futures[i])
    end
    for i = 1:np-nw-1
        int = (i-1)*M_min+nw*(M_min+1)+1+M_min:i*M_min+nw*(M_min+1)+M_min
        u⁺_list[int], u⁻_list[int] = fetch(futures[nw+i])
    end
    
    println("Parallell measurements done. Processing.")
    
    return u⁺_list, u⁻_list
end
