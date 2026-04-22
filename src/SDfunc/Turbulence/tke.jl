export tke_settings, implicit_diffuse!, turbulent_droplet_diffusion!
export calculate_buoyancy_frequency, mixing_length_deardorff!
export duynkerke_tke_timestep!, calculate_Ri, S2_flow_deformation
export no_prog_tke_timestep!, calculate_Pr, smag_lilly_timestep!, smag_lilly_ca_timestep!
export find_bl_height


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
end


function mixing_length_deardorff!(l::Vector{FT}, grid, l_inf::FT, constants) where FT<:AbstractFloat
    c_n = FT(0.76)
    nz  = grid.nz
    for k in 1:nz
        N2   = calculate_buoyancy_frequency(grid, k, constants)
        if N2 > 0 # check sign
            l[k] = min(grid.dz,c_n * sqrt(grid.states.e[k] / N2))
            l[k] = max(l[k], 1e-3)
        else
            l[k] = grid.dz
            l[k] = max(l[k], 1e-3)
        end
    end
    # l[1] = 0.5 * l[2]  # Surface boundary condition
    # l[end] = 0.5 * l[end-1]  # Top boundary condition
end

# function mixing_length_deardorff!(l::Vector{FT}, grid, l_inf::FT, constants) where FT<:AbstractFloat
#     ce1 = 0.19
#     ce2 = 0.51
#     ca = 0.1
#     c1 = 0.76^2 * ca
#     c2 = ce2 + c1
#     nz  = grid.nz
#     for k in 1:nz
#         N2   = calculate_buoyancy_frequency(grid, k, constants)
#         Ri = calculate_Ri(grid, k, constants)
#         if N2 > 0 # check sign
#             l[k] = grid.dz * ( c1 * (1/Ri - 1) - ce1) / c2
#         else
#             l[k] = grid.dz
#         end
#     end
#     # l[1] = 0.5 * l[2]  # Surface boundary condition
#     # l[end] = 0.5 * l[end-1]  # Top boundary condition
# end

function mixing_length_nonlocal!(l::Vector{FT}, grid, l_inf::FT, constants) where FT<:AbstractFloat
    vonkarman = 0.4
    nz  = grid.nz
    for k in 1:nz
        z = grid.centers_z[k]
        len = 30 + 270 * (1 - exp(-z/100)) # example nonlocal length scale
        lc = 1/(vonkarman* z) + 1/len
        l[k] = 1/lc
    end
end

function S2_flow_deformation(grid, k, constants)
    if k == 1 || k == grid.nz
        del_z = grid.dz
    else
        del_z = 2 * grid.dz
    end
    du_dz = (grid.wind.u[min(k+1, grid.nz)] - grid.wind.u[max(k-1, 1)]) / del_z
    dv_dz = (grid.wind.v[min(k+1, grid.nz)] - grid.wind.v[max(k-1, 1)]) / del_z
    S2 = du_dz^2 + dv_dz^2
    return S2
end

function calculate_Ri(grid, k, constants)
    N2 = calculate_buoyancy_frequency(grid, k, constants)
    S2 = S2_flow_deformation(grid, k, constants)
    Ri = N2 / S2
    return Ri
end

function calculate_buoyancy_frequency(grid, k, constants) #where FT<:AbstractFloat
    g = constants.gconst
    qv = grid.states.qv
    ql = grid.states.ql_tmp
    θ = grid.states.θ
    nz = grid.nz
    dz = grid.dz
    # ql = compute_ql_at_cell(grid.states, k)


    kp = min(k + 1, nz)
    km = max(k - 1, 1)
    dz_buo = k == nz ? dz : 2.0 * dz
    dz_buo = k == 1 ? dz : dz_buo
    θv_k = θ[k] * (1 + 0.61 * qv[k] - ql[k])
    θv_up = θ[kp] * (1 + 0.61 * qv[kp] - ql[kp])
    θv_dn = θ[km] * (1 + 0.61 * qv[km] - ql[km])
    N2 = (g / θv_k) * (θv_up - θv_dn) / (dz_buo)

    return N2
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
    d  = zeros(FT, nz)   # main diagonal
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
    max_dθdz = -Inf
    h_idx = nz
    for k in 2:nz
        dθ_dz = (grid.states.θ[k] - grid.states.θ[k-1]) / dz
        if dθ_dz > max_dθdz
            max_dθdz = dθ_dz
            h_idx = k
        end
    end
    h = (h_idx - 1) * dz   # height of inversion base (lower face of cell h_idx)
    return h_idx, h
end

# Non-local counter-gradient correction for θ and qv (Holtslag & Moeng 1991).
# Adds the flux-divergence tendency -∂(K_h γ)/∂z to cells within the boundary layer.
# γ_θ = c_γ * H_kin / (w* h),  w* = (g/θ_sfc * H_v * h)^(1/3)
# Only active when the surface virtual heat flux H_v > 0 (unstable/convective).
function apply_counter_gradient!(grid, tke::tke_settings{FT}, K_h::Vector{FT}, constants, dt::FT) where FT
    (tke.SHF == 0 && tke.LHF == 0) && return

    nz  = grid.nz
    dz  = FT(grid.dz)
    ρ   = grid.states.ρ
    g   = FT(constants.gconst)

    H_kin = tke.SHF / (ρ[1] * constants.Cp_air)          # kinematic sensible heat flux  [K m/s]
    E_kin = tke.LHF / (ρ[1] * constants.L)                # kinematic moisture flux        [kg/kg m/s]
    H_v   = H_kin + FT(0.61) * grid.states.θ[1] * E_kin  # virtual heat flux

    H_v <= 0 && return

    h_idx, h = find_bl_height(grid)
    h <= dz  && return

    w_star = cbrt(g / grid.states.θ[1] * H_v * h)

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

function no_prog_tke_timestep!(grid, tke::tke_settings{FT}, constants, dt::FT) where FT
    nz = grid.nz
    dz = FT(grid.dz)
    e  = grid.states.e
    ρ = grid.states.ρ .+ 0

    l   = zeros(FT, nz)
    K_m = zeros(FT, nz)
    K_h = zeros(FT, nz)
        # grid.area_per_grid
    delta = grid.dz

    grid.states.θ
    grid.states.qv
    # grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, collect(1:nz))

    surface_zonal_momentum_flux = -grid.wind.u[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    surface_meridional_momentum_flux = - grid.wind.v[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)


    # mixing_length_blackadar!(l, grid.centers_z, tke.l_inf)
    mixing_length_deardorff!(l, grid, tke.l_inf, constants)
    # mixing_length_nonlocal!(l, grid, tke.l_inf, constants)

    # Pr = 0.42 #this is dziekan 2019, find else?

    

    for k in 1:nz
        sqrte  = sqrt(max(e[k], tke.e_min))
        K_m[k] = tke.c_m * l[k] * sqrte
        # S2 = S2_flow_deformation(grid, k, constants)
        Pr = calculate_Pr(grid, k, constants)

        # K_h[k] = K_m[k] / Pr 
        K_h[k] = K_m[k] *(1 + 2 * l[k] /delta) #deardorff 1980

    end
    

    #maybe should do rho*qv and rho*theta diffusion to be more consistent with mass conservation just like advection?
    # ρθ = grid.states.ρ .* grid.states.θ
    # ρqv = grid.states.ρ .* grid.states.qv
    # ρu = grid.states.ρ .* grid.wind.u
    # ρv = grid.states.ρ .* grid.wind.v
    implicit_diffuse!(grid.states.θ,  K_h, dt, dz, nz)
    implicit_diffuse!(grid.states.qv, K_h, dt, dz, nz)
    # implicit_diffuse!(ρθ,  K_h, dt, dz, nz)
    # implicit_diffuse!(ρqv, K_h, dt, dz, nz)
    # implicit_diffuse!(grid.states.ρ, K_h, dt, dz, nz)
    implicit_diffuse!(grid.wind.u, K_m, dt, dz, nz,sfc_flux = surface_zonal_momentum_flux)
    implicit_diffuse!(grid.wind.v, K_m, dt, dz, nz, sfc_flux = surface_meridional_momentum_flux)
    # grid.states.θ .= ρθ ./ grid.states.ρ
    # grid.states.qv .= ρqv ./ grid.states.ρ
    # grid.wind.u .= ρu ./ grid.states.ρ
    # grid.wind.v .= ρv ./ grid.states.ρ
    #calculate θl and qtot after diffusion using θ and qv
    grid.states.ρ .= ρ_calc_θ(grid.states.P, grid.states.θ, grid.states.qv, constants)
    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, collect(1:nz))

    apply_counter_gradient!(grid, tke, K_h, constants, dt)

    coriolis_parameter = 2 * constants.Ω * sind(31.5)
    du = zeros(FT, nz)
    dv = zeros(FT, nz)
    for k in 1:nz
        z = grid.centers_z[k]
        du[k] = coriolis_parameter * (grid.wind.v[k] - tke.geostrophic_v(z))
        dv[k] = -coriolis_parameter * (grid.wind.u[k] - tke.geostrophic_u(z))
    end

    grid.wind.u .+= du .* dt
    grid.wind.v .+= dv .* dt


    de = zeros(FT, nz)
    g  = FT(constants.gconst)
    
    P_shear = zeros(FT, nz)
    P_buoy = zeros(FT, nz)
    P_diss = zeros(FT, nz)
    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, collect(1:nz))

    for k in 1:nz
        #Shear
        P_shear[k] = K_m[k] * S2_flow_deformation(grid, k, constants)

        # Buoyancy
        P_buoy[k] = - K_h[k] * calculate_buoyancy_frequency(grid, k, constants)

        # Dissipation
        # c_eps = 1.9*tke.c_m + (0.93 - 1.9*tke.c_m) * l[k] / delta
        c_eps = 0.19 + 0.51 * (l[k]/delta)
        # c_eps = 0.7
        P_diss[k]= c_eps * e[k]^(3/2) / l[k]

        # de[k] = P_shear + P_buoy - ε
        # de[k] =P_buoy
    end
    de .= P_shear + P_buoy - P_diss

    # for k in 1:nz
    #     c_eps      = 0.19 + 0.51 * (l[k] / delta)
    #     diss_coeff = c_eps * sqrt(max(e[k], tke.e_min)) / l[k]
    #     source     = K_m[k] * S2_flow_deformation(grid, k, constants) -
    #                  K_h[k] * calculate_buoyancy_frequency(grid, k, constants)
    #     e[k]       = (e[k] + dt * source) / (1 + dt * diss_coeff)
    # end
    # e .= max.(e, tke.e_min)


    e .+= dt .* de
    e .= max.(e, tke.e_min)

    

    ρe = ρ .* e
    implicit_diffuse!(ρe, K_h, dt, dz, nz)#; sfc_value = e_sfc)
    e .= ρe ./ grid.states.ρ
    # e[1] = max.(e[1], tke.u_star^2 * tke.c_ε)
    e[1] = tke.u_star^2 / sqrt(tke.c_m)
    e .= max.(e, tke.e_min)
    grid.states.e .= e

    return nothing
end



function smag_lilly_timestep!(grid, tke::tke_settings{FT}, constants, dt::FT) where FT
    nz = grid.nz
    dz = FT(grid.dz)
    
    K_m = zeros(FT, nz)
    K_h = zeros(FT, nz)

    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, collect(1:nz))

    surface_zonal_momentum_flux = grid.wind.u[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    surface_meridional_momentum_flux = grid.wind.v[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    
    c_s = 0.17
    delta = grid.dz

    for k in 1:nz
        N2 = calculate_buoyancy_frequency(grid, k, constants)
        S2 = S2_flow_deformation(grid, k, constants)
        Pr = 0.42#calculate_Pr(grid, k, constants)
        fb = N2 > 0 ? max(0, (1 - N2 /(Pr* S2))^1/2) : 1.0
        K_m[k] = c_s^2 * delta^2 * sqrt(S2) * fb
        K_h[k] = K_m[k] / Pr 
    end
    



    implicit_diffuse!(grid.states.θ,  K_h, dt, dz, nz)
    implicit_diffuse!(grid.states.qv, K_h, dt, dz, nz)
    implicit_diffuse!(grid.wind.u, K_h, dt, dz, nz,sfc_flux = surface_zonal_momentum_flux)
    implicit_diffuse!(grid.wind.v, K_h, dt, dz, nz, sfc_flux = surface_meridional_momentum_flux)
    #calculate θl and qtot after diffusion using θ and qv
    grid.states.ρ .= ρ_calc_θ(grid.states.P, grid.states.θ, grid.states.qv, constants)
    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, collect(1:nz))

    coriolis_parameter = 2 * constants.Ω * sind(31.5)
    du = zeros(FT, nz)
    dv = zeros(FT, nz)
    for k in 1:nz
        z = grid.centers_z[k]
        du[k] = coriolis_parameter * (grid.wind.v[k] - tke.geostrophic_v(z))
        dv[k] = -coriolis_parameter * (grid.wind.u[k] - tke.geostrophic_u(z))
    end

    grid.wind.u .+= du .* dt
    grid.wind.v .+= dv .* dt

    return nothing
end


function smag_lilly_ca_timestep!(grid, tke::tke_settings{FT}, constants, dt::FT) where FT
    nz = grid.nz
    dz = FT(grid.dz)

    K_m = zeros(FT, nz)
    K_h = zeros(FT, nz)

    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, collect(1:nz))

    ρ     = grid.states.ρ .+ 0
    speed = sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    surface_zonal_momentum_flux      = ρ[1] * grid.wind.u[1] * tke.u_star^2 / speed
    surface_meridional_momentum_flux = ρ[1] * grid.wind.v[1] * tke.u_star^2 / speed

    c_smag = tke.c_smag
    Pr_t   = tke.Pr_t
    delta  = cbrt(tke.delta_h^2 * dz)   # geometric mean filter scale (Δx*Δy*Δz)^(1/3)

    for k in 1:nz
        N2     = calculate_buoyancy_frequency(grid, k, constants)
        S2     = S2_flow_deformation(grid, k, constants)
        S_norm = sqrt(max(S2, FT(0)))
        Ri     = N2 / (S2 + eps(FT))
        fb     = Ri ≤ 0 ? FT(1) : max(FT(0), 1 - Ri / Pr_t)^(FT(1/4))  # Lilly (1967) stability correction
        νt     = c_smag^2 * (delta * fb)^2 * S_norm
        K_m[k] = νt
        K_h[k] = νt / Pr_t
    end

    # diffuse ρ-weighted quantities: ∂(ρϕ)/∂t = ∂/∂z(ρ K ∂ϕ/∂z)
    ρθ  = ρ .* grid.states.θ
    ρqv = ρ .* grid.states.qv
    ρu  = ρ .* grid.wind.u
    ρv  = ρ .* grid.wind.v

    implicit_diffuse!(ρθ,  K_h, dt, dz, nz)
    implicit_diffuse!(ρqv, K_h, dt, dz, nz)
    implicit_diffuse!(ρu,  K_m, dt, dz, nz, sfc_flux = surface_zonal_momentum_flux)
    implicit_diffuse!(ρv,  K_m, dt, dz, nz, sfc_flux = surface_meridional_momentum_flux)

    grid.states.θ .= ρθ  ./ ρ
    grid.states.qv .= ρqv ./ ρ
    grid.wind.u   .= ρu  ./ ρ
    grid.wind.v   .= ρv  ./ ρ

    grid.states.ρ .= ρ_calc_θ(grid.states.P, grid.states.θ, grid.states.qv, constants)
    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, collect(1:nz))

    coriolis_parameter = 2 * constants.Ω * sind(31.5)
    du = zeros(FT, nz)
    dv = zeros(FT, nz)
    for k in 1:nz
        z = grid.centers_z[k]
        du[k] =  coriolis_parameter * (grid.wind.v[k] - tke.geostrophic_v(z))
        dv[k] = -coriolis_parameter * (grid.wind.u[k] - tke.geostrophic_u(z))
    end
    grid.wind.u .+= du .* dt
    grid.wind.v .+= dv .* dt

    return nothing
end


function turbulent_droplet_diffusion!(droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT, coagdata, constants) where FT
    #edited from grabowski and abade 2017 eqn 10, Dziekan 2019 uses this too.
    #their version builds on w' from the previous timestep
    nz    = grid.nz
    dz    = FT(grid.dz)
    Z_max = FT(nz) * dz
    Ct = 1.5
    l  = zeros(FT, nz)
    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, 1:nz)
    # θ = theta_from_T(grid.states.T, grid.states.P, constants)
    θ = grid.states.θ
    mixing_length_deardorff!(l, grid, tke.l_inf,constants)
    # mixing_length_nonlocal!(l, grid, tke.l_inf, constants)

    for z in 1:nz
        if droplets.grid_range[z] == nothing
            continue
        end
        
        k = coagdata.I[droplets.grid_range[z]]
        tau = l[z]/((2pi)^(1/3)) * (Ct/grid.states.e[z])^(1/2) 
        sigma = sqrt(grid.states.e[z] * 2/3)

        w_prime = sqrt(1- exp(-2*dt/tau)) * sigma
        droplets.z_loc[k] .+= w_prime .* rand(Uniform(-0.1,0.1),length(k)) * dt
    end
    
    

    return
end
