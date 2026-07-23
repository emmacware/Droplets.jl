export S2_flow_deformation#, calculate_Ri, calculate_buoyancy_frequency,
export calculate_buoyancy_frequency
export implicit_diffuse!
#moist_buoyancy_production,
# calculate_Pr, find_bl_height, apply_counter_gradient!
export tke_settings#, AbstractTurbulenceScheme, ProgTKE, MellorYamada, SmagorinskyLilly, KNMIturb
# export AbstractMixingLengthFunction, Deardorff, SL, Nonlocal, Blackadar, KNMI
export turbulent_droplet_diffusion!, turbulent_droplet_diffusion_exact!
export dθldz, moist_buoyancy_frequency

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

struct tke_settings{FT<:AbstractFloat, GU, GV}
    e_min::FT
    u_star::FT
    SHF::FT
    LHF::FT
    geostrophic_u::GU
    geostrophic_v::GV
    dry_buoyancy::Bool
    c_n::FT
    vk::FT
    bott_α::FT
    bott_β::FT
    my_diss::FT
    GH_lims::Tuple{FT,FT}
    c_m::FT  # EDMF-X TKE eddy-diffusivity coefficient, Lopez-Gomez et al. (2020) Table 1
    c_d::FT  # EDMF-X TKE dissipation coefficient, Lopez-Gomez et al. (2020) Table 1
    c_b::FT  # EDMF-X static-stability mixing-length coefficient (CLIMAParameters default)
end

function tke_settings{FT}(;
        e_min::FT               = FT(1e-3),
        u_star::FT              = FT(0.25),
        SHF::FT                 = FT(16),
        LHF::FT                 = FT(93),
        geostrophic_u           = z -> FT(0),
        geostrophic_v           = z -> FT(0),
        dry_buoyancy::Bool      = false,
        c_n::FT                 = FT(0.76),
        vk::FT                  = FT(0.4),
        bott_α::FT              = FT(0.05),
        bott_β::FT              = FT(4.0),
        my_diss::FT             = FT(1/16.6),
        GH_lims::Tuple{FT,FT}   = (FT(-0.4), FT(0.0233)),
        c_m::FT                 = FT(0.14),
        c_d::FT                 = FT(0.22),
        c_b::FT                 = FT(0.4),
    ) where {FT<:AbstractFloat}
    tke_settings{FT, typeof(geostrophic_u), typeof(geostrophic_v)}(
        e_min, u_star, SHF, LHF, geostrophic_u, geostrophic_v,
        dry_buoyancy, c_n, vk, bott_α, bott_β, my_diss, GH_lims, c_m, c_d, c_b)
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

    RdCp = Rd / Cp
    LCp  = L / Cp
    Π_k  = (P[k]  / P0)^RdCp;  θlv_k  = (θ[k]  - LCp * ql[k]  / Π_k)  * (1 + 0.61 * qv[k]  - ql[k])
    Π_kp = (P[kp] / P0)^RdCp;  θlv_kp = (θ[kp] - LCp * ql[kp] / Π_kp) * (1 + 0.61 * qv[kp] - ql[kp])
    Π_km = (P[km] / P0)^RdCp;  θlv_km = (θ[km] - LCp * ql[km] / Π_km) * (1 + 0.61 * qv[km] - ql[km])
    N2 = (g / θlv_k) * (θlv_kp - θlv_km) / dz_buo
    return N2
end

function moist_buoyancy_frequency(grid, k, constants)
    g  = constants.gconst
    nz = grid.nz
    dz = grid.dz
    kp = min(k + 1, nz)
    km = max(k - 1, 1)
    dz_buo = (k == 1 || k == nz) ? dz : 2.0 * dz

    L  = constants.L
    Cp = constants.Cp_air
    Cp_vapor = constants.Cp_vapor
    Cp_water = constants.Cp_water
    Rd = constants.Rd
    P0 = constants.P0
    qv = grid.states.qv
    ql = grid.states.ql_tmp
    θ  = grid.states.θ
    T = grid.states.T_tmp
    P  = grid.states.P

    # ------------------------------------------------------------
    # Total water (NO singular transform)
    # ------------------------------------------------------------
    qt  = qv[k] + ql[k]
    qt_p = qv[kp] + ql[kp]
    qt_m = qv[km] + ql[km]

    # ------------------------------------------------------------
    # Virtual temperature approximation (Cotton/DK consistent)
    # ------------------------------------------------------------
    Tv = T[k] * (1 + 0.61 * qv[k] - ql[k])

    # ------------------------------------------------------------
    # Virtual potential temperature
    # ------------------------------------------------------------
    θv = θ[k] * (1 + 0.61 * qv[k] - ql[k])

    # ------------------------------------------------------------
    # Moist adiabatic lapse rate correction (DK-style form)
    # ------------------------------------------------------------
    es = esat(T[k])
    rsat = constants.ϵ * es / (P[k] - es)

    # moist static stability correction term (cleaned form)
    Γd = g / Cp
    denom = 1 + (L * rsat) / (Rd * T[k])

    Γm = Γd / denom

    # ------------------------------------------------------------
    # θv gradient (core DK structure)
    # ------------------------------------------------------------
    dθv_dz = (θv - θ[k] * (1 + 0.61 * qv[k] - ql[k])) / dz_buo

    # fallback centered form (more consistent numerically)
    dθv_dz = (θ[kp]*(1 + 0.61*qv[kp] - ql[kp]) -
              θ[km]*(1 + 0.61*qv[km] - ql[km])) / dz_buo

    # ------------------------------------------------------------
    # Moisture gradient correction (Cotton-style)
    # ------------------------------------------------------------
    dqt_dz = (qt_p - qt_m) / dz_buo

    # DK-style buoyancy frequency
    N2 = (g / θv) * dθv_dz - (g / (1 + qt)) * dqt_dz

    # lapse rate correction (optional DK enhancement)
    N2 += Γm * (g / T[k])
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
    dt::FT, dz::FT, nz::Int,turbdata::turbulence_data;
    sfc_flux::Union{FT,Nothing} = nothing) where FT

    #zero out the arrays in turbulence_data
    turbdata.K_faces .= 0
    turbdata.rhs .= 0
    turbdata.a .= 0
    turbdata.b .= 0
    turbdata.c .= 0
    turbdata.d .= 0

    # Face diffusivities
    K_face = turbdata.K_faces#zeros(FT, nz + 1)

    for k in 2:nz
        K_face[k] = FT(0.5) * (K_centers[k-1] + K_centers[k])
    end
    K_face[1] = 0 # (bottom NoFlux)
    K_face[nz+1] = 0 # (top NoFlux)

    r = dt / (dz^2)
    rhs = turbdata.rhs#zeros(FT, nz)   # right-hand side (copy of ϕ)
    rhs .= ϕ

    a = turbdata.a#zeros(FT, nz)   # sub-diagonal
    b = turbdata.b#zeros(FT, nz)   # main diagonal
    c = turbdata.c#zeros(FT, nz)   # super-diagonal
    d = turbdata.d#zeros(FT, nz)   # RHS

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

function turbulent_droplet_diffusion!(::DynOFF,l,droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
end
# function turbulent_droplet_diffusion!(::DynON,l,droplets::droplet_attributes_1d{FT},
#     grid, tke::tke_settings{FT}, dt::FT) where FT
#     #edited from grabowski and abade 2018 to be just dependent on TKE
#     #also using gillespie solution for dw'
#     nz    = grid.nz
#     dz    = grid.dz

#     Ct = FT(0.63)
#     Ce = FT(0.89)
#     n_sub  = 8  # substep the OU integration within this outer dt; e/l/tau stay fixed across substeps
#     dt_sub = dt / n_sub
#     inv_idx = findfirst(k -> grid.states.qv[k] < 0.008, 1:grid.nz)
#     z_inv   = isnothing(inv_idx) ? FT(Inf) : FT((inv_idx - 1) * grid.dz)

#     for z in 1:nz
#         isempty(droplets.grid_range[z]) && continue
#         k_all = droplets.I[droplets.grid_range[z]]

#         r = volume_to_radius.(droplets.X[k_all])
#         cloud_mask = r .< FT(40e-6)   # only cloud droplets < 40 μm
#         k = k_all[cloud_mask]
#         isempty(k) && continue

#         zm = max(z - 1, 1); zp = min(z + 1, nz)
#         e_eff = (grid.states.e[zm] + grid.states.e[z] + grid.states.e[zp]) / 3
#         l_z   = (l[zm] + l[z] + l[zp]) / 3

#         tau = Ct*l_z * sqrt(Ce / e_eff)

#         sigma2    = e_eff * FT(2/3)
#         noise_amp = sqrt(2*sigma2 / tau)
#         w_corr    = 1 - dt_sub / (2*tau)
#         dt_sqrt   = sqrt(dt_sub)
#         dt_32     = dt_sub * dt_sqrt

#         w0     = copy(droplets.w_prime[k])
#         z_disp = zeros(FT, length(k))

#         for _ in 1:n_sub
#             n1 = randn(FT, length(k))
#             n2 = randn(FT, length(k))

#             dz_disp = w0 .* (dt_sub * w_corr) .+ noise_amp .* FT(0.5) .* (n1 .+ n2 ./ sqrt(FT(3))) .* dt_32
#             dw      = (-w0 .* (dt_sub / tau) .+ noise_amp .* n1 .* dt_sqrt) .* w_corr

#             z_disp .+= dz_disp
#             w0       = w0 .+ dw
#         end

#         droplets.z_loc[k]   .+= z_disp
#         droplets.z_loc[k]   .= min.(droplets.z_loc[k], z_inv)  # capped at the inversion: no crossing up
#         droplets.w_prime[k]  .= w0

#         droplets.cell_id[k] = clamp.(floor.(Int, droplets.z_loc[k] / grid.dz) .+ 1, -1, nz)
#     end

    
#     return
# end

# function turbulent_droplet_diffusion!(::DynON,K_h,droplets::droplet_attributes_1d{FT},
#     grid, tke::tke_settings{FT}, dt::FT) where FT

#     nz    = grid.nz
    
#     for z in 1:nz
#         isempty(droplets.grid_range[z]) && continue
#         k     = droplets.I[droplets.grid_range[z]]
#         dz_d = z == 1 ? grid.dz : z == nz ? grid.dz : 2*grid.dz
#         dKdz = (K_h[min(z+1, nz)] - K_h[max(z-1, 1)]) / dz_d
#         Δz = dKdz.*dt .+ sqrt(2*K_h[z]*dt).*randn(FT, length(k))
#         droplets.z_loc[k] .+= min(Δz,5)
#         droplets.cell_id[k] = clamp.(floor.(Int,droplets.z_loc[k]/ grid.dz) .+ 1, -1, nz)
        
#     end
#     return
#     sort!(droplets.I, by = i -> droplets.cell_id[i])
# end

function turbulent_droplet_diffusion!(::DynON,l,droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
    nz    = grid.nz
    dz    = grid.dz
    Ct = FT(0.63)
    Ce = FT(0.89)
    e = grid.states.e
    for z in 1:nz
        isempty(droplets.grid_range[z]) && continue
        k     = droplets.I[droplets.grid_range[z]]
        e_z   = max(e[z], tke.e_min)
        l_z   = dz
        tau   = Ct * l_z * sqrt(Ce / e_z)
        sigma2 = e_z * FT(2/3)
        # Thomson (1987) well-mixed correction: ½ ∂σ_w²/∂z prevents spurious
        # accumulation in low-TKE regions (e.g. above the inversion)
        e_lo  = max(z > 1  ? e[z-1] : e[z], tke.e_min)
        e_hi  = max(z < nz ? e[z+1] : e[z], tke.e_min)
        dsigma2_dz = FT(2/3) * (e_hi - e_lo) / (2 * dz)
        w_amp = sqrt((1 - exp(-2 * dt / tau)) * sigma2)
        droplets.w_prime[k] .*= exp(-dt / tau)
        droplets.w_prime[k] .+= FT(0.5) * dsigma2_dz * dt
        droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))
    end
    return
end

# function turbulent_droplet_diffusion!(::DynON,l,droplets::droplet_attributes_1d{FT},
#     grid, tke::tke_settings{FT}, dt::FT) where FT
#     #edited from grabowski and abade 2018 eqn 10, 
#     nz    = grid.nz
#     # l_z = grid.dz
#     Ct = FT(0.63)
#     C_e = FT(0.89)
#     for z in 1:nz
#         l_z = min(l[z], grid.dz)
#         isempty(droplets.grid_range[z]) && continue
#         k     = droplets.I[droplets.grid_range[z]]
#         e_z   = grid.states.e[z]
#         # eps_z = tke.my_diss * (2 * max(e_z, tke.e_min))^(FT(3/2)) / l_z
#         # eps_z = C_e * max(e_z, tke.e_min)^(FT(3/2)) / l_z

#         tau   = Ct * l_z * (C_e / max(e_z, tke.e_min))^(FT(1/2))

#         sigma = sqrt(e_z * FT(2/3))
#         w_amp = sqrt(2 * sigma^2/tau)
#         dw = -droplets.w_prime[k]./tau .* dt 
#         droplets.w_prime[k] .+= dw
#         droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))* sqrt(dt)

#     end

# Weil et al. (2004) Eq. (15) Lagrangian stochastic model for SGS droplet transport.
#
# SDE:  dw' = [-w'/τ_t + ½(1/σ²_w · dσ²_w/dt · w' + ∂σ²_w/∂z)] dt + √(2σ²_w/τ_t · dt) η
#
# where:
#   σ²_w = (2/3) e              (vertical velocity variance, isotropic TKE)
#   ε    = B₁⁻¹ (2e)^(3/2) / l (MY dissipation rate)
#   τ_t  = 2σ²_w / (C_o ε)     (Lagrangian timescale, C_o = 3)
#
# Drift terms (Thomson 1987 well-mixed condition):
#   temporal: ½ (dσ²_w/dt)/σ²_w · w' ≈ ½ (de/dt)/e · w'  — pass e_prev to include
#   spatial:  ½ ∂σ²_w/∂z = (1/3) ∂e/∂z                   — always included
#
# Discretisation: exact OU for the decay+noise, Euler for the drift corrections.
function weil_turbulent_droplet_diffusion!(::DynOFF, l, droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT, e_prev=nothing) where FT
end

function weil_turbulent_droplet_diffusion!(::DynON, l, droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT, e_prev=nothing) where FT

    nz     = grid.nz
    dz     = FT(grid.dz)
    Co     = FT(3)        # Kolmogorov constant C_o (Weil et al. 2004)
    n_sub  = 10
    dt_sub = dt / n_sub
    e      = grid.states.e

    for z in 1:nz
        isempty(droplets.grid_range[z]) && continue
        k = droplets.I[droplets.grid_range[z]]

        e_z = max(e[z], tke.e_min)
        l_z = l[z]

        # σ²_w = (2/3)e,   ε = B₁⁻¹(2e)^(3/2)/l
        sigma2 = FT(2/3) * e_z
        eps_z  = max(tke.my_diss * (2*e_z)^FT(1.5) / l_z, FT(1e-10))

        # Lagrangian timescale τ_t = 2σ²_w / (C_o ε)
        tau = max(2*sigma2 / (Co * eps_z), dt_sub)

        # Spatial drift: ½ ∂σ²_w/∂z = (1/3) ∂e/∂z  [m/s per substep]
        kp = min(z+1, nz);  km = max(z-1, 1)
        dz_eff     = (z == 1 || z == nz) ? dz : 2*dz
        dsigma2_dz = FT(2/3) * (e[kp] - e[km]) / dz_eff
        drift_z    = FT(0.5) * dsigma2_dz * dt_sub

        # Temporal drift coefficient: ½ (de/dt)/e  [s⁻¹, multiplied into w_prime per substep]
        temporal_coeff = if isnothing(e_prev)
            FT(0)
        else
            FT(0.5) * (e_z - max(e_prev[z], tke.e_min)) / (e_z * dt)
        end

        # Pre-compute per-substep constants (τ, σ², gradients fixed across substeps)
        decay = exp(-dt_sub / tau)
        w_amp = sqrt(sigma2 * (1 - exp(-2*dt_sub / tau)))  # exact OU noise amplitude

        Z_max = nz * dz
        for _ in 1:n_sub
            droplets.w_prime[k] .*= (decay + temporal_coeff * dt_sub)  # OU decay + temporal drift
            droplets.w_prime[k] .+= drift_z                             # spatial drift
            droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))
            droplets.z_loc[k]   .+= droplets.w_prime[k] .* dt_sub
            # Reflecting boundary conditions: flip velocity and reflect z_loc
            for ki in k
                if droplets.z_loc[ki] > Z_max
                    droplets.z_loc[ki]   = 2*Z_max - droplets.z_loc[ki]
                    droplets.w_prime[ki] = -droplets.w_prime[ki]
                elseif droplets.z_loc[ki] < 0
                    droplets.z_loc[ki]   = -droplets.z_loc[ki]
                    droplets.w_prime[ki] = -droplets.w_prime[ki]
                end
            end
            droplets.cell_id[k] = clamp.(floor.(Int, droplets.z_loc[k] / dz) .+ 1, 1, nz)
        end
    end
    return
end

#     return
# end