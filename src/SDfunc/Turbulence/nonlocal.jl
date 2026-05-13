export S2_flow_deformation#, calculate_Ri, calculate_buoyancy_frequency,
export calculate_buoyancy_frequency
#moist_buoyancy_production,
# implicit_diffuse!,# calculate_Pr, find_bl_height, apply_counter_gradient!
export tke_settings#, AbstractTurbulenceScheme, ProgTKE, MellorYamada, SmagorinskyLilly, KNMIturb
# export AbstractMixingLengthFunction, Deardorff, SL, Nonlocal, Blackadar, KNMI
export turbulent_droplet_diffusion!, turbulent_droplet_diffusion_exact!
export dθldz

# abstract type AbstractTurbulenceScheme end
# struct ProgTKE <: AbstractTurbulenceScheme end
# struct MellorYamada <: AbstractTurbulenceScheme end
# struct SmagorinskyLilly <: AbstractTurbulenceScheme end
# struct KNMIturb <: AbstractTurbulenceScheme end

# abstract type AbstractMixingLengthFunction end
# struct Deardorff <: AbstractMixingLengthFunction end
# struct SL <: AbstractMixingLengthFunction end
# struct Nonlocal <: AbstractMixingLengthFunction end
# struct Blackadar <: AbstractMixingLengthFunction end
# struct KNMI <: AbstractMixingLengthFunction end

Base.@kwdef struct tke_settings{FT<:AbstractFloat}
    e_min::FT   = FT(1e-12)     # TKE floor (m²/s²)
    u_star::FT  = FT(0.25)     # surface friction velocity (m/s)
    SHF::FT     = FT(16)      # surface sensible heat flux (W/m²)
    LHF::FT     = FT(93)      # surface latent heat flux (W/m²)
    geostrophic_u::Function = z -> FT(0)  # geostrophic u wind profile (m/s)
    geostrophic_v::Function = z -> FT(0)  # geostrophic v wind profile (m/s)
    # turbulence_scheme::AbstractTurbulenceScheme = ProgTKE()
    # mixing_length_function::AbstractMixingLengthFunction = Deardorff()
    dry_buoyancy::Bool = false           # use dry θ in N²; false = moist θlv
    c_n ::FT = FT(0.76)                     # Deardorff mixing length constant
    vk::FT = FT(0.4)                        # von Karman constant
    bott_α::FT = FT(0.05)                     # Bott mixing length parameter
    bott_β::FT = FT(4.0)                     # Bott mixing length parameter
    my_diss::FT = FT(1/16.6)                       # constant TKE dissipation
    GH_lims::Tuple{FT,FT} = (FT(-0.4), FT(0.03))  
end


function S2_flow_deformation(grid, k, constants)
    if k == 1
        du_dz = 0.5*(grid.wind.u[2] - grid.wind.u[1]) / grid.dz
        dv_dz = 0.5*(grid.wind.v[2] - grid.wind.v[1]) / grid.dz
        # du_dz = (grid.wind.u[1]) / 0.5*grid.dz
        # dv_dz = (grid.wind.v[1]) / 0.5*grid.dz
    elseif k == grid.nz
        du_dz = (grid.wind.u[end] - grid.wind.u[end-1]) / grid.dz
        dv_dz = (grid.wind.v[end] - grid.wind.v[end-1]) / grid.dz
    else
        du_dz = (grid.wind.u[k+1] - grid.wind.u[k-1]) / (2*grid.dz)
        dv_dz = (grid.wind.v[k+1] - grid.wind.v[k-1]) / (2*grid.dz)
    end
    S2 = du_dz^2 + dv_dz^2
    return S2
end

# function calculate_Ri(grid, k, constants)
#     N2 = calculate_buoyancy_frequency(grid, k, constants)
#     S2 = S2_flow_deformation(grid, k, constants)
#     Ri = N2 / S2
#     return Ri
# end

function calculate_buoyancy_frequency(grid, k, constants; dry::Bool=false)
    g  = constants.gconst
    nz = grid.nz
    dz = grid.dz
    kp = min(k + 1, nz)
    km = max(k - 1, 1)
    dz_buo = (k == 1 || k == nz) ? dz : 2.0 * dz

    if dry
        θ = grid.states.θ
        return (g / θ[k]) * (θ[kp] - θ[km]) / dz_buo
    end

    L  = constants.L
    Cp = constants.Cp_air
    Rd = constants.Rd
    P0 = constants.P0
    qv = grid.states.qv
    ql = grid.states.ql_tmp
    θ  = grid.states.θ
    P  = grid.states.P

    function θlv(kk)
        Π  = (P[kk] / P0)^(Rd / Cp)
        θl = θ[kk] - (L / Cp) * ql[kk] / Π
        return θl * (1 + 0.608 * qv[kk] - ql[kk])
    end

    θlv_k  = θlv(k)
    N2 = (g / θlv_k) * (θlv(kp) - θlv(km)) / dz_buo
    return N2
end

function dθldz(θl_t, grid, k, constants)
    nz = grid.nz
    dz = grid.dz
    kp = min(k + 1, nz)
    km = max(k - 1, 1)
    dz_buo = (k == 1 || k == nz) ? dz : 2.0 * dz


    return (θl_t[kp] - θl_t[km]) / dz_buo

end





function implicit_diffuse!(ϕ::Vector{FT}, K_centers::Vector{FT},
    dt::FT, dz::FT, nz::Int;
    sfc_flux::Union{FT,Nothing} = nothing) where FT

    # Face diffusivities
    K_face = zeros(FT, nz + 1)
    for k in 2:nz
        K_face[k] = FT(0.5) * (K_centers[k-1] + K_centers[k])
    end
    K_face[1] = 0 # (bottom NoFlux)
    K_face[nz+1] = 0 # (top NoFlux)

    r = dt / (dz^2)
    rhs = ϕ .+ 0 

    a = zeros(FT, nz)   # sub-diagonal
    b = zeros(FT, nz)   # main diagonal
    c = zeros(FT, nz)   # super-diagonal
    d = zeros(FT, nz)   # RHS

    for k in 2:nz-1
        Km = K_face[k]       
        Kp = K_face[k+1]
        a[k] = -r/2 * Km
        b[k] = 1 + r/2 * (Km + Kp)
        c[k] = -r/2 * Kp
        d[k] = ϕ[k] + r/2 * (Kp*(ϕ[k+1]-ϕ[k]) - Km*(ϕ[k]-ϕ[k-1]))
    end
    
    # b[1] = 1; d[1] = ϕ[1]
    # Tridiagonal coefficients for k=1 (first interior cell)
    a[1] = 0               # no sub-diagonal for first cell
    b[1] = 1 + r/2*(K_face[1]+K_face[2])
    c[1] = -r/2 * K_face[2]
    d[1] = ϕ[1] + r/2 * (K_face[2]*(ϕ[2]-ϕ[1]))#- r/2 * K_face[1]*(sfc_flux * dt / dz) 
    a[nz] = -r/2 * K_face[nz]
    b[nz] = 1 + r/2 * K_face[nz]
    c[nz] = 0
    d[nz] = ϕ[nz] - r/2 * K_face[nz] * (ϕ[nz] - ϕ[nz-1])
    if sfc_flux !== nothing
        d[1] += dt * (sfc_flux / (dz))
    end


    c_star = zeros(FT, nz)
    d_star = zeros(FT, nz)
    c_star[1] = c[1] / b[1]
    d_star[1] = d[1] / b[1]
    bet = b[1]

    for k in 2:nz
        # c_star[k] = c[k-1]/bet
        
        bet = b[k] - a[k] * c_star[k-1]
        c_star[k] = c[k] / bet
        d_star[k] = (d[k] - a[k] * d_star[k-1]) / bet
    end


    rhs[nz] = d_star[nz]
    for k in nz-1:-1:1
        rhs[k] = d_star[k] - c_star[k] * rhs[k+1]
    end

    ϕ .= rhs
    

    return nothing
end

# function Fc(Ri)
#     if Ri < 0
#         return sqrt(1 - 18 * Ri )
#     else
#         return 1 / (1 + 10 * Ri * ( 1 + 8 * Ri))
#     end
# end


# function calculate_Pr(grid, k, constants) #from Nishizawa 2015
#     Ri = calculate_Ri(grid, k, constants)
#     if Ri < 0
#         Pr= 0.7 * sqrt((1 - 16 *Ri)/(1-40*Ri))
#     elseif 0< Ri < 0.25
#         Pr= 0.7 * 1/(1-(1-0.7)*Ri/0.25)
#     else
#         Pr = 1
#     end
#     Pr = min(Pr, 2.0)
#     return Pr
# end


# function find_bl_height(grid)
#     nz = grid.nz
#     dz = grid.dz
#     h_idx = nz
#     #inversion is where qv drops below 8 g/kg
#     for k in 1:nz
#         if grid.states.qv[k] < 0.008
#             h_idx = k
#             break
#         end
#     end
#     h = grid.centers_z[h_idx] 
#     # height of inversion base (lower face of cell h_idx)
#     return h_idx, h
# end

# Non-local counter-gradient correction for θ and qv (Holtslag & Moeng 1991).
# Adds the flux-divergence tendency -∂(K_h γ)/∂z to cells within the boundary layer.
# γ_θ = c_γ * H_kin / (w* h),  w* = (g/θ_sfc * H_v * h)^(1/3)
function apply_counter_gradient!(grid, tke::tke_settings{FT}, K_h::Vector{FT}, constants, dt::FT) where FT
    (tke.SHF == 0 && tke.LHF == 0) && return

    nz  = grid.nz
    dz  = FT(grid.dz)
    ρ   = grid.states.ρ
    g   = FT(constants.gconst)

    HF = tke.SHF / (ρ[1] * constants.Cp_air)        
    LHF = tke.LHF / (ρ[1] * constants.L)                
    H_v   = HF + FT(0.61) * grid.states.θ[1] * LHF  

    H_v <= 0 && return

    h_idx, h = find_bl_height(grid)
    h <= dz  && return
    
    #buoyancy parameter Beta_g
    Beta_g = g / grid.states.θ[1]

    w_star = cbrt(g / grid.states.θ[1] * HF * h)

    γ_θ  = tke.c_γ * H_kin / (w_star * h)
    γ_qv = tke.c_γ * E_kin / (w_star * h)

    # Tendency: ∂ϕ/∂t|_nl = -∂(K_h γ)/∂z
    # Non-local flux K_h γ is zero at the surface (k=1 lower face) and at the BL top (h_idx upper face).
    for k in 1:h_idx-1
        K_up = k < h_idx-1 ? FT(0.5)*(K_h[k] + K_h[k+1]) : FT(0)
        K_dn = k > 1       ? FT(0.5)*(K_h[k-1] + K_h[k]) : FT(0)
        grid.states.θ[k]  -= dt * (K_up - K_dn) * γ_θ  / dz
        grid.states.qv[k] -= dt * (K_up - K_dn) * γ_qv / dz
    end
end


function turbulent_droplet_diffusion!(l,droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
    #edited from grabowski and abade 2017 eqn 10, Dziekan 2019 uses this too.
    nz    = grid.nz
    
    Ct = FT(0.63)
    for z in 1:nz
        droplets.grid_range[z] === nothing && continue
        k     = droplets.I[droplets.grid_range[z]]
        e_z   = max(grid.states.e[z],0)
        l_z = min(grid.dz, l[z])
        tau   = l_z / (2π)^(FT(1/3)) * sqrt(Ct / e_z)
        # @assert dt < 2.0 * tau
        # tau = grid.dz / sqrt(e_z) # eddy turnover time based on grid scale and local TKE
        sigma = sqrt(e_z * FT(2/3))
        w_amp = sqrt(1 - exp(-2 * dt / tau)) * sigma
        droplets.w_prime[k] .*= exp(-dt / tau)
        droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))
    end

    return
end


# function turbulent_droplet_diffusion!(l,droplets::droplet_attributes_1d{FT},
#     grid, tke::tke_settings{FT}, dt::FT) where FT
#     #edited from grabowski and abade 2018 eqn 10, 
#     nz    = grid.nz
#     # l_z = grid.dz
#     Ct = FT(0.63)
#     C_e = FT(0.89)
#     for z in 1:nz
#         l_z = min(l[z], grid.dz)
#         droplets.grid_range[z] === nothing && continue
#         k     = droplets.I[droplets.grid_range[z]]
#         e_z   = grid.states.e[z]
#         eps_z = tke.my_diss * (2 * max(e_z, tke.e_min))^(FT(3/2)) / l_z
#         # eps_z = C_e * max(e_z, tke.e_min)^(FT(3/2)) / l_z

#         tau   = Ct * (l_z )^(2/3) / eps_z^(1/3)

#         sigma = sqrt(e_z * FT(2/3))
#         w_amp = sqrt(2 * sigma^2/tau)
#         dw = -droplets.w_prime[k]./tau .* dt 
#         droplets.w_prime[k] .+= dw
#         droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))* sqrt(dt)

#     end

#     return
# end


# function turbulent_droplet_diffusion!(l, droplets::droplet_attributes_1d{FT},
#     grid, tke::tke_settings{FT}, dt::FT) where FT
#     # OU process with Thomson (1987) well-mixed correction.
#     # Without the ∂σ²/∂z term, droplets accumulate in low-TKE regions (spurious downward drift).
#     nz  = grid.nz
#     dz  = FT(grid.dz)
#     Ct  = FT(0.63)
#     ce  = tke.my_diss
#     e   = grid.states.e
#     for z in 1:nz
#         droplets.grid_range[z] === nothing && continue
#         k     = droplets.I[droplets.grid_range[z]]
#         l_z   = l[z]
#         e_z   = max(e[z], tke.e_min)
#         eps_z = ce * (2 * e_z)^(FT(3/2)) / l_z
#         tau   = Ct * l_z^(FT(2/3)) / eps_z^(FT(1/3))
#         sigma2  = e_z * FT(2/3)
#         decay   = exp(-dt / tau)
#         w_amp   = sqrt(sigma2 * (1 - decay^2))

#         # e_up = z < nz ? max(e[z+1], tke.e_min) : e_z
#         # e_dn = z >  1 ? max(e[z-1], tke.e_min) : e_z
#         # dz_eff         = (z == 1 || z == nz) ? dz : 2*dz
#         # drift_correction = FT(0.5) * FT(2/3) * (e_up - e_dn) / dz_eff * dt

#         droplets.w_prime[k] .*= decay
#         # droplets.w_prime[k] .+= drift_correction
#         droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))
#     end
#     return
# end