export tke_settings, tke_timestep!, implicit_diffuse!, turbulent_droplet_diffusion!


Base.@kwdef struct tke_settings{FT<:AbstractFloat}
    c_m::FT    = FT(0.1)      # momentum diffusivity constant
    c_ε::FT    = FT(0.93)     # dissipation constant
    l_inf::FT  = FT(70.0)     # Blackadar asymptotic mixing length (m)
    e_min::FT  = FT(1e-6)     # TKE floor (m²/s²)
    u_star::FT = FT(0.25)     # surface friction velocity (m/s)
    geostrophic_u::Function = z -> FT(0)  # geostrophic u wind profile (m/s)
    geostrophic_v::Function = z -> FT(0)  # geostrophic v wind profile (m/s)
end

# function mixing_length_blackadar!(l::Vector{FT}, z_centers, l_inf::FT,
#                                   κ::FT = FT(0.4)) where FT
#     for k in eachindex(z_centers)
#         z    = max(FT(z_centers[k]), FT(0.5))   # avoid z = 0
#         l[k] = κ * z / (1 + κ * z / l_inf)
#     end
# end


# function mixing_length_deardorff!(l::Vector{FT}, grid, l_inf::FT, z_surf::FT = FT(0)) where FT
#     c_n = 0.76
#     #delta is the geometric mean of the grid spacing, Dziekan 2019 just uses dz
#     delta = grid.dz
#     for k in eachindex(grid.centers_z)
#         N2 = calculate_buoyancy_frequency(grid,k)
#         if N2 > 0
#             l[k] = min(delta, c_n * sqrt(grid.states.e[k]/N2))
#         else
#             l[k] = delta
#         end
#     end
# end
function mixing_length_deardorff!(l::Vector{FT}, grid, l_inf::FT, z_surf::FT = FT(0)) where FT
    c_n = FT(0.76)
    κ   = FT(0.4)
    for k in eachindex(grid.centers_z)
        z    = max(FT(grid.centers_z[k]) - z_surf, FT(0.5))
        l_bl = κ * z / (1 + κ * z / l_inf)   # Blackadar asymptotic length
        N2   = calculate_buoyancy_frequency(grid, k)
        if N2 > 0
            l[k] = min(l_bl, c_n * sqrt(grid.states.e[k] / N2))
        else
            l[k] = l_bl
        end
    end
end

function calculate_buoyancy_frequency(grid, k)
    g = constants.gconst
    Tv_k = T_virtual(grid.states.T[k], grid.states.qv[k])
    Tv_kp = T_virtual(grid.states.T[min(k+1, grid.nz)], grid.states.qv[min(k+1, grid.nz)])
    θv_k = theta_from_T(Tv_k, grid.states.P[k], constants)
    θv_kp = theta_from_T(Tv_kp, grid.states.P[min(k+1, grid.nz)], constants)

    N2 = g * (θv_kp - θv_k) / (θv_k * grid.dz)
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
    # K_face[1] = 0  (bottom NoFlux)
    # K_face[nz+1] = 0  (top NoFlux)

    r = dt / (dz^2)
    d  = zeros(FT, nz)   # main diagonal
    rhs = copy(ϕ)
    a = -r/2
    b = 1 + r
    c = -r/2


    a = zeros(FT, nz)   # sub-diagonal
    b = zeros(FT, nz)   # main diagonal
    c = zeros(FT, nz)   # super-diagonal
    d = zeros(FT, nz)   # RHS

    # Interior rows
    for k in 1:nz
        Km = K_face[k]       
        Kp = K_face[k+1]
        a[k] = -r/2 * Km
        b[k] = 1 + r/2 * (Km + Kp)
        c[k] = -r/2 * Kp
        d[k] = ϕ[k] + r/2 * (Kp*(ϕ[k+1 > nz ? k : k+1]-ϕ[k]) - Km*(ϕ[k]-ϕ[k-1 < 1 ? k : k-1]))
    end

    if sfc_value !== nothing
        d[1] = sfc_value
    end

    b[1] = 1; d[1] = ϕ[1]
    b[nz] = 1; d[nz] = ϕ[nz]

    c_star = zeros(FT, nz-1)
    d_star = zeros(FT, nz)
    c_star[1] = c[1] / b[1]
    d_star[1] = d[1] / b[1]
    for k in 2:nz-1
        denom = b[k] - a[k] * c_star[k-1]
        c_star[k] = c[k] / denom
    end

    for k in 2:nz
        numerator = d[k] - a[k] * d_star[k-1]
        denom = b[k] - a[k] * c_star[k-1]
        d_star[k] = numerator / denom
    end

    rhs[nz] = d_star[nz]
    for k in nz-1:-1:1
        rhs[k] = d_star[k] - c_star[k] * rhs[k+1]
    end
    rhs[end] = ϕ[end]
    rhs[1] = ϕ[1]
    ϕ .= rhs
    

    if sfc_flux !== nothing
        ϕ[1] += dt * (sfc_flux / dz)
    end


end





function tke_timestep!(grid, tke::tke_settings{FT}, constants, dt::FT) where FT
    nz = grid.nz
    dz = FT(grid.dz)
    e  = grid.states.e

    l   = zeros(FT, nz)
    K_m = zeros(FT, nz)
    K_h = zeros(FT, nz)
        # grid.area_per_grid
    delta = grid.dz


    # mixing_length_blackadar!(l, grid.centers_z, tke.l_inf)
    mixing_length_deardorff!(l, grid, tke.l_inf)

    for k in 1:nz
        sqrte  = sqrt(max(e[k], tke.e_min))
        K_m[k] = tke.c_m * l[k] * sqrte
        # K_h[k] = (1 + 2*l[k])/delta * K_m[k]
        K_h[k] = K_m[k] / 0.42 # Prandtl number of 0.42, Dziekan 2019

    end
    

    coriolis_parameter = 2 * constants.Ω * sin(deg2rad(31.5))
    du = zeros(FT, nz)
    dv = zeros(FT, nz)
    for k in 1:nz
        z = grid.centers_z[k]
        du[k] = coriolis_parameter * (grid.wind.v[k] - tke.geostrophic_v(z))
        dv[k] = -coriolis_parameter * (grid.wind.u[k] - tke.geostrophic_u(z))
    end

    surface_zonal_momentum_flux = grid.wind.u[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    surface_meridional_momentum_flux = grid.wind.v[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)

    implicit_diffuse!(grid.states.θ,  K_h, dt, dz, nz)
    implicit_diffuse!(grid.states.qv, K_h, dt, dz, nz)
    implicit_diffuse!(grid.wind.u, K_h, dt, dz, nz, sfc_flux = surface_zonal_momentum_flux)
    implicit_diffuse!(grid.wind.v, K_h, dt, dz, nz, sfc_flux = surface_meridional_momentum_flux)


    de = zeros(FT, nz)
    g  = FT(constants.gconst)

    for k in 2:nz
        #Shear
        du_dz = (grid.wind.u[k] - grid.wind.u[k-1]) / dz
        dv_dz = (grid.wind.v[k] - grid.wind.v[k-1]) / dz
        P_shear = K_m[k] * (du_dz^2 + dv_dz^2)

        # Buoyancy
        kp = min(k + 1, nz)
        dz_buo = k == nz ? dz : FT(2) * dz
        Tv_up = T_virtual(grid.states.T[kp],  grid.states.qv[kp])
        Tv_dn = T_virtual(grid.states.T[k-1], grid.states.qv[k-1])
        Tv_k = T_virtual(grid.states.T[k],   grid.states.qv[k])
        θv_k = theta_from_T(Tv_k, grid.states.P[k], constants)
        θv_up = theta_from_T(Tv_up, grid.states.P[kp], constants)
        θv_dn = theta_from_T(Tv_dn, grid.states.P[k-1], constants)
        P_buoy = -(g / θv_k) * K_h[k] * (θv_up - θv_dn) / dz_buo

        # Dissipation
        c_eps = 1.9*tke.c_m + (0.93 - 1.9*tke.c_m) * l[k] / delta
        ε= c_eps * max(e[k], FT(0))^FT(1.5) / l[k]

        de[k] = P_shear + P_buoy - ε
    end

    e[2:end] .+= dt .* de[2:end]
    e .= max.(e, tke.e_min)
    e_sfc = FT(tke.u_star^2 / sqrt(tke.c_m)) # Ackerman 1995 says this plus 0.35*w_star^2, 0 under stable conditions
    #transport
    implicit_diffuse!(e, K_m, dt, dz, nz)#; sfc_value = e_sfc)
    e .= max.(e, tke.e_min)

    grid.states.T .= T_from_theta(grid.states.θ, grid.states.P, constants)
    grid.states.ρ .= ρ_ideal_gas(grid.states.P, grid.states.T, grid.states.qv, constants)

    return nothing
end


function turbulent_droplet_diffusion!(droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT,coagdata) where FT
    #edited from grabowski and abade 2017 eqn 10, Dziekan 2019 uses this too.
    #their version builds on w' from the previous timestep 
    nz    = grid.nz
    dz    = FT(grid.dz)
    Z_max = FT(nz) * dz
    Ct = 1.5
    l  = zeros(FT, nz)
    mixing_length_deardorff!(l, grid, tke.l_inf)

    for z in 1:nz
        if droplets.grid_range[z] == nothing
            continue
        end
        
        k = coagdata.I[droplets.grid_range[z]]
        tau = l[z]/((2pi)^(1/3)) * (Ct/grid.states.e[z])^(1/2) 
        sigma = sqrt(grid.states.e[z] * 2/3)

        w_prime = sqrt(1- exp(-2*dt/tau)) * sigma
        droplets.z_loc[k] .+= w_prime .* randn(length(k)) * dt
    end
    
    

    return
end
