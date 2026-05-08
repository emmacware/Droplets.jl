export S2_flow_deformation, calculate_Ri, calculate_buoyancy_frequency,
moist_buoyancy_production,
implicit_diffuse!, implicit_diffuse_theta!, calculate_Pr, find_bl_height, apply_counter_gradient!
export tke_settings, AbstractTurbulenceScheme, ProgTKE, MellorYamada, SmagorinskyLilly, KNMIturb
export AbstractMixingLengthFunction, Deardorff, SL, Nonlocal, Blackadar, KNMI

abstract type AbstractTurbulenceScheme end
struct ProgTKE <: AbstractTurbulenceScheme end
struct MellorYamada <: AbstractTurbulenceScheme end
struct SmagorinskyLilly <: AbstractTurbulenceScheme end
struct KNMIturb <: AbstractTurbulenceScheme end

abstract type AbstractMixingLengthFunction end
struct Deardorff <: AbstractMixingLengthFunction end
struct SL <: AbstractMixingLengthFunction end
struct Nonlocal <: AbstractMixingLengthFunction end
struct Blackadar <: AbstractMixingLengthFunction end
struct KNMI <: AbstractMixingLengthFunction end

Base.@kwdef struct tke_settings{FT<:AbstractFloat}
    c_m::FT     = FT(0.1)      # momentum diffusivity constant
    c_ε::FT     = FT(0.93)     # dissipation constant
    l_inf::FT   = FT(70.0)     # Blackadar asymptotic mixing length (m)
    e_min::FT   = FT(1e-9)     # TKE floor (m²/s²)
    u_star::FT  = FT(0.25)     # surface friction velocity (m/s)
    c_smag::FT  = FT(0.17)     # Smagorinsky constant
    Pr_t::FT    = FT(1/3)      # turbulent Prandtl number (neutral stratification)
    delta_h::FT = FT(100.0)    # horizontal grid spacing for filter scale (m)
    SHF::FT     = FT(16)      # surface sensible heat flux (W/m²)
    LHF::FT     = FT(93)      # surface latent heat flux (W/m²)
    c_γ::FT     = FT(7.2)      # Holtslag & Moeng (1991) counter-gradient constant
    geostrophic_u::Function = z -> FT(0)  # geostrophic u wind profile (m/s)
    geostrophic_v::Function = z -> FT(0)  # geostrophic v wind profile (m/s)
    turbulence_scheme::AbstractTurbulenceScheme = ProgTKE()
    mixing_length_function::AbstractMixingLengthFunction = Deardorff()
    A_e::FT = FT(0.2)                    # entrainment efficiency (Moeng 1986)
    dry_buoyancy::Bool = false           # use dry θ in N²; false = moist θlv
end


function S2_flow_deformation(grid, k, constants)
    if k == 1
        du_dz = (grid.wind.u[2] - grid.wind.u[1]) / grid.dz
        dv_dz = (grid.wind.v[2] - grid.wind.v[1]) / grid.dz
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

function calculate_Ri(grid, k, constants)
    N2 = calculate_buoyancy_frequency(grid, k, constants)
    S2 = S2_flow_deformation(grid, k, constants)
    Ri = N2 / S2
    return Ri
end

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

# Buoyancy production term using L&H (2004) eq. 13 moist A,B decomposition.
# Unsaturated: A = 1+0.61qt, B = 0.61
# Saturated:   A and B from moist adiabatic correction (accounts for latent heat release).
function moist_buoyancy_production(grid, k::Int, K_h_k, constants)
    g  = constants.gconst
    L  = constants.L
    Cp = constants.Cp_air
    Rd = constants.Rd
    P0 = constants.P0
    nz = grid.nz
    dz = grid.dz

    function cell_vars(kk)
        P  = grid.states.P[kk]
        θ  = grid.states.θ[kk]
        qv = grid.states.qv[kk]
        ql = grid.states.ql_tmp[kk]
        Π  = (P / P0)^(Rd / Cp)
        θl = θ - (L / Cp) * ql / Π
        qt = qv + ql
        θlv = θl * (1 + 0.608 * qv - ql)
        return θl, qt, θlv, θ * Π  # T = θ * Π
    end

    kp = min(k + 1, nz)
    km = max(k - 1, 1)
    dz_eff = (k == 1 || k == nz) ? dz : 2dz

    θl_up, qt_up, _,      _   = cell_vars(kp)
    θl_dn, qt_dn, _,      _   = cell_vars(km)
    θl_k,  qt_k,  θlv_k, T_k = cell_vars(k)

    dθl_dz = (θl_up - θl_dn) / dz_eff
    dqt_dz = (qt_up - qt_dn) / dz_eff

    ql_k = grid.states.ql_tmp[k]
    if ql_k > 0
        qs   = grid.states.qv[k]   # qv = qs at saturation
        fac  = 0.622 * L / (Rd * T_k)
        facp = L / (T_k * Cp)
        A = (1 - qt_k + 1.61 * qs * (1 + fac)) / (1 + fac * facp * qs)
        B = facp * A - 1
    else
        A = 1 + 0.61 * qt_k
        B = 0.61
    end

    N2_moist = (g / θlv_k) * (A * dθl_dz + B * θl_k * dqt_dz)
    return -K_h_k * N2_moist
end


function implicit_diffuse!(ϕ::Vector{FT}, K_centers::Vector{FT},
    dt::FT, dz::FT, nz::Int;
    sfc_value::Union{FT,Nothing} = nothing,
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

    # rhs[end] = ϕ[end]
    # rhs[1] = ϕ[1]
    ϕ .= rhs#rhs
    

    # if sfc_flux !== nothing
    #     ϕ[1] += dt * (sfc_flux / dz)
    # end

    return nothing
end

# Generalized Crank-Nicolson diffusion solver (Paegle et al. 1976).
# β_cur + β_fut = 1; paper uses β_cur=0.75, β_fut=0.25 (θ-method with θ=β_fut).
# Standard CN: β_cur=β_fut=0.5. Fully implicit: β_cur=0, β_fut=1.
# WARNING: β_fut < 0.5 is only conditionally stable (need K*dt/dz² < 1/(1-2*β_fut)).
function implicit_diffuse_theta!(ϕ::Vector{FT}, K_centers::Vector{FT},
    dt::FT, dz::FT, nz::Int;
    β_cur::FT = FT(0.5),
    β_fut::FT = FT(0.5),
    sfc_flux::Union{FT,Nothing} = nothing) where FT

    K_face = zeros(FT, nz + 1)
    for k in 2:nz
        K_face[k] = FT(0.5) * (K_centers[k-1] + K_centers[k])
    end
    K_face[1]    = FT(0)
    K_face[nz+1] = FT(0)

    r = dt / dz^2

    a = zeros(FT, nz)
    b = zeros(FT, nz)
    c = zeros(FT, nz)
    d = zeros(FT, nz)

    # Bottom cell (NoFlux below: K_face[1]=0)
    Kp = K_face[2]
    a[1] = FT(0)
    b[1] = 1 + β_fut * r * Kp
    c[1] = -β_fut * r * Kp
    d[1] = ϕ[1] + β_cur * r * Kp * (ϕ[2] - ϕ[1])

    # Interior cells
    for k in 2:nz-1
        Km = K_face[k]
        Kp = K_face[k+1]
        a[k] = -β_fut * r * Km
        b[k] =  1 + β_fut * r * (Km + Kp)
        c[k] = -β_fut * r * Kp
        d[k] =  ϕ[k] + β_cur * r * (Kp * (ϕ[k+1] - ϕ[k]) - Km * (ϕ[k] - ϕ[k-1]))
    end

    # Top cell (NoFlux above: K_face[nz+1]=0)
    Km = K_face[nz]
    a[nz] = -β_fut * r * Km
    b[nz] =  1 + β_fut * r * Km
    c[nz] = FT(0)
    d[nz] =  ϕ[nz] - β_cur * r * Km * (ϕ[nz] - ϕ[nz-1])

    if sfc_flux !== nothing
        d[1] += dt * sfc_flux / dz
    end

    # Thomas algorithm (forward sweep)
    c_star = zeros(FT, nz)
    d_star = zeros(FT, nz)
    c_star[1] = c[1] / b[1]
    d_star[1] = d[1] / b[1]
    for k in 2:nz
        bet      = b[k] - a[k] * c_star[k-1]
        c_star[k] = c[k] / bet
        d_star[k] = (d[k] - a[k] * d_star[k-1]) / bet
    end

    # Back substitution
    ϕ[nz] = d_star[nz]
    for k in nz-1:-1:1
        ϕ[k] = d_star[k] - c_star[k] * ϕ[k+1]
    end

    return nothing
end

function Fc(Ri)
    if Ri < 0
        return sqrt(1 - 18 * Ri )
    else
        return 1 / (1 + 10 * Ri * ( 1 + 8 * Ri))
    end
end


function calculate_Pr(grid, k, constants) #from Nishizawa 2015
    Ri = calculate_Ri(grid, k, constants)
    if Ri < 0
        Pr= 0.7 * sqrt((1 - 16 *Ri)/(1-40*Ri))
    elseif 0< Ri < 0.25
        Pr= 0.7 * 1/(1-(1-0.7)*Ri/0.25)
    else
        Pr = 1
    end
    Pr = min(Pr, 2.0)
    return Pr
end


function find_bl_height(grid)
    nz = grid.nz
    dz = grid.dz
    h_idx = nz
    #inversion is where qv drops below 8 g/kg
    for k in 1:nz
        if grid.states.qv[k] < 0.008
            h_idx = k
            break
        end
    end
    h = grid.centers_z[h_idx] 
    # height of inversion base (lower face of cell h_idx)
    return h_idx, h
end

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