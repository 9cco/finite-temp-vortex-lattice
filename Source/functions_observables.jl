####################################################################################################################
#                            Planar Structure function
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
# Returns the un-normalized local vorticity by preforming a plaquette sum using the gauge-invariant
# difference of the θ field.
function n⁺(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
	return (mod(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺, two_pi) - ϕ.A[1] + mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₁.θ⁺, two_pi) - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos)  
			- mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₂.θ⁺, two_pi) + ϕᵣ₊₂.A[1]
			- mod(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺, two_pi) + (ϕ.A[2] + two_pi*c.f*(h_pos-1)))
end
function n⁻(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
    return (mod(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻, two_pi) - ϕ.A[1] + mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₁.θ⁻, two_pi) - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos) 
        - mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₂.θ⁻, two_pi) + ϕᵣ₊₂.A[1] 
        - mod(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻, two_pi) + (ϕ.A[2] + two_pi*c.f*(h_pos-1)))
end
# Returns both vorticities of the gauge-invariant phase difference.
function nᵣGaugeInv(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
    vort_θ⁺ = (mod(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺ - ϕ.A[1], two_pi) + mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₁.θ⁺ - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos), two_pi)
        - mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₂.θ⁺ - ϕᵣ₊₂.A[1], two_pi) 
        - mod(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺ - (ϕ.A[2] + two_pi*c.f*(h_pos-1)), two_pi))
    vort_θ⁻ = (mod(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻ - ϕ.A[1], two_pi) + mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₁.θ⁻ - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos), two_pi)
        - mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₂.θ⁻ - ϕᵣ₊₂.A[1], two_pi) 
        - mod(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻ - (ϕ.A[2] + two_pi*c.f*(h_pos-1)), two_pi))
    return vort_θ⁺, vort_θ⁻
end
#function nᵣ(c::SystConstants, ϕ::LatticeSite, ϕᵣ₊₁::LatticeSite, ϕᵣ₊₂::LatticeSite, ϕᵣ₊₁₊₂::LatticeSite, h_pos::Int64)
#    vort_θ⁺ = (mod(ϕᵣ₊₁.θ⁺ - ϕ.θ⁺, two_pi) - ϕ.A[1] + mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₁.θ⁺, two_pi) - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos)
#        - mod(ϕᵣ₊₁₊₂.θ⁺ - ϕᵣ₊₂.θ⁺, two_pi)  + ϕᵣ₊₂.A[1]
#        - mod(ϕᵣ₊₂.θ⁺ - ϕ.θ⁺, two_pi)  + (ϕ.A[2] + two_pi*c.f*(h_pos-1)))
#    vort_θ⁻ = (mod(ϕᵣ₊₁.θ⁻ - ϕ.θ⁻, two_pi) - ϕ.A[1] + mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₁.θ⁻, two_pi) - (ϕᵣ₊₁.A[2] + two_pi*c.f*h_pos)
#        - mod(ϕᵣ₊₁₊₂.θ⁻ - ϕᵣ₊₂.θ⁻, two_pi)  + ϕᵣ₊₂.A[1]
#        - mod(ϕᵣ₊₂.θ⁻ - ϕ.θ⁻, two_pi)  + (ϕ.A[2] + two_pi*c.f*(h_pos-1)))
#    return vort_θ⁺, vort_θ⁻
#end
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


# -----------------------------------------------------------------------------------------------------------
# Assuming we have a state ψ, we want to find the lattice of vortexes.
function vortexSnapshot(ψ::State)
    L = ψ.consts.L
    V⁺ = zeros(L,L)
    V⁻ = zeros(L,L)
    
    # Sum over the corners
    # Upper left corner
     ϕ = ψ.lattice[1,1]
     ϕᵣ₊₁ = ψ.lattice[1,2]
     ϕᵣ₊₂ = ψ.lattice[L,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,2]
    (V⁺[1,1], V⁻[1,1]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)
    
    # Lower left corner
     ϕ = ψ.lattice[L,1]
     ϕᵣ₊₁ = ψ.lattice[L,2]
     ϕᵣ₊₂ = ψ.lattice[L-1,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,2]
    (V⁺[L,1], V⁻[L,1]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)
    
    # Lower right corner
     ϕ = ψ.lattice[L,L]
     ϕᵣ₊₁ = ψ.lattice[L,1]
     ϕᵣ₊₂ = ψ.lattice[L-1,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,1]
    (V⁺[L,L], V⁻[L,L]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)
    
    # Upper right corner
     ϕ = ψ.lattice[1,L]
     ϕᵣ₊₁ = ψ.lattice[1,1]
     ϕᵣ₊₂ = ψ.lattice[L,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,1]
    (V⁺[1,L], V⁻[1,L]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)
    
    # Sum over borders without corners
    for i = 2:(L-1)
        # Upper border (counting in +x direction)
        ϕ = ψ.lattice[1,i]
        ϕᵣ₊₁ = ψ.lattice[1,i+1]
        ϕᵣ₊₂ = ψ.lattice[L,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L,i+1]
        (V⁺[1,i], V⁻[1,i]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)
        
        # Left border (counting in -y direction)
        ϕ = ψ.lattice[i,1]
        ϕᵣ₊₁ = ψ.lattice[i,2]
        ϕᵣ₊₂ = ψ.lattice[i-1,1]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,2]
        (V⁺[i,1], V⁻[i,1]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)
        
        # Lower border (counting in +x direction)
        ϕ = ψ.lattice[L,i]
        ϕᵣ₊₁ = ψ.lattice[L,i+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L-1,i+1]
        (V⁺[L,i], V⁻[L,i]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)
        
        # Right border (counting in -y direction)
        ϕ = ψ.lattice[i,L]
        ϕᵣ₊₁ = ψ.lattice[i,1]
        ϕᵣ₊₂ = ψ.lattice[i-1,L]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,1]
        (V⁺[i,L], V⁻[i,L]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)
    end
    
    # Sum over rest of the bulk
    for h_pos = 2:(L-1)
        for v_pos = 2:(L-1)
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos,h_pos+1]
            ϕᵣ₊₂ = ψ.lattice[v_pos-1,h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[v_pos-1,h_pos+1]
            (V⁺[v_pos,h_pos], V⁻[v_pos, h_pos]) = nᵣ(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)
        end
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

# -----------------------------------------------------------------------------------------------------------
# V⁺ and V⁻ are LxL matrices of plaquettes containing the vorticities nᵣ of the two components.
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
    
    # Sum over the corners
    # Upper left corner
     r = [0, L-1] # For r we assume origo is in position [L,1] of the lattice. 
                  # Note that r is the same as pos (found previously) with y-axis flipped and -1 in each direction.
                  # Additionally we define it such that we get the usual r = [x,y] order of dimensions.
     ϕ = ψ.lattice[1,1]
     ϕᵣ₊₁ = ψ.lattice[1,2]
     ϕᵣ₊₂ = ψ.lattice[L,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,2]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower left corner
     r = [0,0]
     ϕ = ψ.lattice[L,1]
     ϕᵣ₊₁ = ψ.lattice[L,2]
     ϕᵣ₊₂ = ψ.lattice[L-1,1]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,2]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
    
    # Lower right corner
     r = [L-1,0]
     ϕ = ψ.lattice[L,L]
     ϕᵣ₊₁ = ψ.lattice[L,1]
     ϕᵣ₊₂ = ψ.lattice[L-1,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L-1,1]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Upper right corner
     r = [L-1,L-1]
     ϕ = ψ.lattice[1,L]
     ϕᵣ₊₁ = ψ.lattice[1,1]
     ϕᵣ₊₂ = ψ.lattice[L,L]
     ϕᵣ₊₁₊₂ = ψ.lattice[L,1]
    sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    
    # Sum over borders without corners
    for i = 2:(L-1)
        # Upper border (counting in +x direction)
        r = [i-1, L-1]
        ϕ = ψ.lattice[1,i]
        ϕᵣ₊₁ = ψ.lattice[1,i+1]
        ϕᵣ₊₂ = ψ.lattice[L,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L,i+1]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Left border (counting in -y direction)
        r = [0, L-i]
        ϕ = ψ.lattice[i,1]
        ϕᵣ₊₁ = ψ.lattice[i,2]
        ϕᵣ₊₂ = ψ.lattice[i-1,1]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,2]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, 1)*exp(im*(k⋅r))
        
        # Lower border (counting in +x direction)
        r = [i-1, 0]
        ϕ = ψ.lattice[L,i]
        ϕᵣ₊₁ = ψ.lattice[L,i+1]
        ϕᵣ₊₂ = ψ.lattice[L-1,i]
        ϕᵣ₊₁₊₂ = ψ.lattice[L-1,i+1]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, i)*exp(im*(k⋅r))
        
        # Right border (counting in -y direction)
        r = [L-1, L-i]
        ϕ = ψ.lattice[i,L]
        ϕᵣ₊₁ = ψ.lattice[i,1]
        ϕᵣ₊₂ = ψ.lattice[i-1,L]
        ϕᵣ₊₁₊₂ = ψ.lattice[i-1,1]
        sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
        sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, L)*exp(im*(k⋅r))
    end
    
    # Sum over rest of the bulk
    for h_pos = 2:(L-1)
        for v_pos = 2:(L-1)
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos,h_pos+1]
            ϕᵣ₊₂ = ψ.lattice[v_pos-1,h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[v_pos-1,h_pos+1]
            sum⁺ += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
            sum⁻ += n⁻(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    
    return (abs2(sum⁺),abs2(sum⁻))
end

# -----------------------------------------------------------------------------------------------------------
# This version should do the same as above, but determines the neighbors dynamically at a cost to performance.
function structureFunctionPlussDyn{T<:Real}(k::Array{T,1}, ψ::State)
    L = ψ.consts.L
    sum = Complex(0.0)
    
    # Sum over entire lattice and determine nearest neighbors dynamically.
    for h_pos = 1:L
        for v_pos = 1:L
            
            r = [h_pos-1, L-v_pos]
            ϕ = ψ.lattice[v_pos,h_pos]
            ϕᵣ₊₁ = ψ.lattice[v_pos, mod(h_pos, L)+1]
            ϕᵣ₊₂ = ψ.lattice[mod(v_pos-2,L)+1, h_pos]
            ϕᵣ₊₁₊₂ = ψ.lattice[mod(v_pos-2,L)+1, mod(h_pos, L)+1]
            sum += n⁺(ψ.consts, ϕ, ϕᵣ₊₁, ϕᵣ₊₂, ϕᵣ₊₁₊₂, h_pos)*exp(im*(k⋅r))
        end
    end
    return abs2(sum)
end


####################################################################################################################
#                            Thermal averages
#
####################################################################################################################

# -----------------------------------------------------------------------------------------------------------
function structureFunctionPlussAvg{T<:Real}(k::Array{T,1}, syst::SystConstants, M::Int64, Δt::Int64)
    S = zeros(M)
    secondMoment = 0
    s_norm = 1/(sim.f*sim.L^2*two_pi)^2
    # Finding a state that has reached thermal equilibrium
    (t₀, dE, ψ, ψ₂) = findEquilibrium(syst)
    
    # Loop over M measurements
    for m = 1:M
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ)
        end
        
        # Make a normalized measurement
        secondMoment = (S[m] = s_norm*structureFunctionPluss(k, ψ))^2
    end
    
    # Return average and error estimates
    av = mean(S)
    secondMoment = secondMoment/M
    τₒ = autocorrTime(S, 5.0)
    err = (1+2*τₒ)*(secondMoment - av^2)/(M-1)
    
    return (av, err, S)
end
# Take in a matrix of k-values and calculate both the vorticity of θ⁺ and θ⁻.
function structureFunctionAvg{T<:Real}(ks::Array{Array{T, 1}, 2}, syst::SystConstants, M::Int64, Δt::Int64)
    Lky = size(ks, 1)
    Lkx = size(ks, 2)
    S⁺ = [zeros(M) for y=1:Lky, x=1:Lkx]    # Matrix containing the series of measurements for each k
    S⁻ = [zeros(M) for y=1:Lky, x=1:Lkx]
    Sm⁺ = [0.0 for y=1:Lky, x=1:Lkx]
    Sm⁻ = [0.0 for y=1:Lky, x=1:Lkx]
    s_norm_inv = 1/(syst.f*syst.L^2*two_pi)^2
    
    # Finding a state that has eached thermal equilibrium
    println("Finding equilibrium")
    (t₀, dE, ψ, ψ₂) = findEquilibrium(syst)
    
    prinln("Making measurements")
    # Loop over M measurements
    for m = 1:M
        print("Measurement progress: $(Int(round(m/M*100,0)))% \r")
        flush(STDOUT)
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ)
        end
        
        # Make measurements 
        for y = 1:Lky, x = 1:Lkx
            (S⁺[y,x][m], S⁻[y,x][m]) = structureFunction(ks[y,x], ψ)
            Sm⁺[y,x] += S⁺[y,x][m]^2
            Sm⁻[y,x] += S⁻[y,x][m]^2
        end
    end
    
    # Return average and error estimates
    avS⁺ = [mean(S⁺[y,x]) for y=1:Lky, x=1:Lkx]
    avS⁻ = [mean(S⁻[y,x]) for y=1:Lky, x=1:Lkx]
    τ⁺ = [autocorrTime(S⁺[y,x], 5.0) for y=1:Lky, x=1:Lkx]
    τ⁻ = [autocorrTime(S⁻[y,x], 5.0) for y=1:Lky, x=1:Lkx]
    errS⁺ = [0.0 for y=1:Lky, x=1:Lkx]
    errS⁻ = [0.0 for y=1:Lky, x=1:Lkx]
    
    for y=1:Lky, x=1:Lkx
        Sm⁺[y,x] /= M
        Sm⁻[y,x] /= M
        errS⁺[y,x] = (1+2*τ⁺[y,x])*(Sm⁺[y,x] - avS⁺[y,x]^2)/(M-1)
        errS⁻[y,x] = (1+2*τ⁻[y,x])*(Sm⁻[y,x] - avS⁻[y,x]^2)/(M-1)
    end
    
    return (avS⁺, errS⁺, avS⁻, errS⁻)
end

# --------------------------------------------------------------------------------------------------
# Take in a matrix of k-values and calculate both the vorticity of θ⁺ and θ⁻.
# Same as previous structureFunction, but now assumes the state is at equilibrium
function structureFunctionAvg!{T<:Real}(ks::Array{Array{T, 1}, 2}, ψ::State, sim::Controls, M::Int64, Δt::Int64)
    syst = ψ.consts
    Lky = size(ks, 1)
    Lkx = size(ks, 2)
    S⁺ = [zeros(M) for y=1:Lky, x=1:Lkx]    # Matrix containing the series of measurements for each k
    S⁻ = [zeros(M) for y=1:Lky, x=1:Lkx]
    Sm⁺ = [0.0 for y=1:Lky, x=1:Lkx]
    Sm⁻ = [0.0 for y=1:Lky, x=1:Lkx]
    s_norm_inv = 1/(syst.f*syst.L^2*two_pi)^2
    
    println("Making measurements over a $(Lkx)×$(Lky) matrix of ks.")
    # Loop over M measurements
    for m = 1:M
        print("Measurement progress: $(Int(round(m/M*100,0)))% \r")
        flush(STDOUT)
        
        # Take Δt MCS
        for i = 1:Δt
            mcSweep!(ψ, sim)
        end
        
        # Make measurements 
        for y = 1:Lky, x = 1:Lkx
            (S⁺[y,x][m], S⁻[y,x][m]) = s_norm_inv.*structureFunction(ks[y,x], ψ)
            Sm⁺[y,x] += S⁺[y,x][m]^2
            Sm⁻[y,x] += S⁻[y,x][m]^2
        end
    end
    
    # Return average and error estimates
    avS⁺ = [mean(S⁺[y,x]) for y=1:Lky, x=1:Lkx]
    avS⁻ = [mean(S⁻[y,x]) for y=1:Lky, x=1:Lkx]
    τ⁺ = [autocorrTime(S⁺[y,x], 5.0) for y=1:Lky, x=1:Lkx]
    τ⁻ = [autocorrTime(S⁻[y,x], 5.0) for y=1:Lky, x=1:Lkx]
    errS⁺ = [0.0 for y=1:Lky, x=1:Lkx]
    errS⁻ = [0.0 for y=1:Lky, x=1:Lkx]
    
    for y=1:Lky, x=1:Lkx
        Sm⁺[y,x] /= M
        Sm⁻[y,x] /= M
        errS⁺[y,x] = (1+2*τ⁺[y,x])*(Sm⁺[y,x] - avS⁺[y,x]^2)/(M-1)
        errS⁻[y,x] = (1+2*τ⁻[y,x])*(Sm⁻[y,x] - avS⁻[y,x]^2)/(M-1)
    end
    
    return (avS⁺, errS⁺, S⁺, avS⁻, errS⁻, S⁻)
end

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

