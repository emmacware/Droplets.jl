export tke_settings, tke_timestep!, implicit_diffuse!, turbulent_droplet_diffusion!


Base.@kwdef struct tke_settings{FT<:AbstractFloat}
    c_m::FT    = FT(0.1)      # momentum diffusivity constant
    c_h::FT    = FT(0.135)    # scalar diffusivity constant  (c_m / Pr_t, Pr_t ≈ 0.74)
    c_ε::FT    = FT(0.93)     # dissipation constant
    l_inf::FT  = FT(70.0)     # Blackadar asymptotic mixing length (m)
    e_min::FT  = FT(1e-6)     # TKE floor (m²/s²)
    u_star::FT = FT(0.25)     # surface friction velocity (m/s)
end

function mixing_length_blackadar!(l::Vector{FT}, z_centers, l_inf::FT,
                                  κ::FT = FT(0.4)) where FT
    for k in eachindex(z_centers)
        z    = max(FT(z_centers[k]), FT(0.5))   # avoid z = 0
        l[k] = κ * z / (1 + κ * z / l_inf)
    end
end

function implicit_diffuse!(ϕ::Vector{FT}, K_centers::Vector{FT},
                           dt::FT, dz::FT, nz::Int;
                           sfc_value::Union{FT,Nothing} = nothing) where FT

    # Face diffusivities
    K_face = zeros(FT, nz + 1)
    for k in 2:nz
        K_face[k] = FT(0.5) * (K_centers[k-1] + K_centers[k])
    end
    # K_face[1] = 0  (bottom NoFlux)
    # K_face[nz+1] = 0  (top NoFlux)

    r = dt / dz^2
    dl  = zeros(FT, nz)   # sub-diagonal
    d   = zeros(FT, nz)   # main diagonal
    du  = zeros(FT, nz)   # super-diagonal
    rhs = copy(ϕ)

    for k in 1:nz
        dl[k] = -r * K_face[k]
        du[k] = -r * K_face[k + 1]
        d[k]  = FT(1) - dl[k] - du[k]
    end

    if sfc_value !== nothing
        dl[1]  = FT(0)
        du[1]  = FT(0)
        d[1]   = FT(1)
        rhs[1] = sfc_value
    end

    for k in 2:nz
        m      = dl[k] / d[k-1]
        d[k]   -= m * du[k-1]
        rhs[k] -= m * rhs[k-1]
    end
    ϕ[nz] = rhs[nz] / d[nz]
    for k in nz-1:-1:1
        ϕ[k] = (rhs[k] - du[k] * ϕ[k+1]) / d[k]
    end
end


function tke_timestep!(grid, tke::tke_settings{FT}, constants, dt::FT) where FT
    nz = grid.nz
    dz = FT(grid.dz)
    e  = grid.states.e

    l   = zeros(FT, nz)
    K_m = zeros(FT, nz)
    K_h = zeros(FT, nz)

    mixing_length_blackadar!(l, grid.centers_z, tke.l_inf)


    for k in 1:nz
        sqrte  = sqrt(max(e[k], tke.e_min))
        K_m[k] = tke.c_m * l[k] * sqrte
        K_h[k] = tke.c_h * l[k] * sqrte
    end


    implicit_diffuse!(grid.states.θ,  K_h, dt, dz, nz)
    implicit_diffuse!(grid.states.qv, K_h, dt, dz, nz)
    # dont diffuse u and v, scm prescribes them


    de = zeros(FT, nz)
    g  = FT(constants.gconst)

    for k in 2:nz
        du_dz = (grid.wind.u[k] - grid.wind.u[k-1]) / dz
        dv_dz = (grid.wind.v[k] - grid.wind.v[k-1]) / dz
        P_shear = K_m[k] * (du_dz^2 + dv_dz^2)

        kp = min(k + 1, nz)
        dz_buo = k == nz ? dz : FT(2) * dz
        Tv_up = T_virtual(grid.states.T[kp],  grid.states.qv[kp])
        Tv_dn = T_virtual(grid.states.T[k-1], grid.states.qv[k-1])
        Tv_k = T_virtual(grid.states.T[k],   grid.states.qv[k])
        P_buoy = -(g / Tv_k) * K_h[k] * (Tv_up - Tv_dn) / dz_buo

        ε= tke.c_ε * max(e[k], FT(0))^FT(1.5) / l[k]

        de[k]   = P_shear + P_buoy - ε
    end

    e[2:end] .+= dt .* de[2:end]
    e .= max.(e, tke.e_min)
    e_sfc = FT(tke.u_star^2 / sqrt(tke.c_m))
    implicit_diffuse!(e, K_m, dt, dz, nz; sfc_value = e_sfc)
    e .= max.(e, tke.e_min)

    grid.states.T .= T_from_theta(grid.states.θ, grid.states.P, constants)
    grid.states.ρ .= ρ_ideal_gas(grid.states.P, grid.states.T, grid.states.qv, constants)

    return nothing
end


function turbulent_droplet_diffusion!(droplets::droplet_attributes_1d{FT},
                                      grid, tke::tke_settings{FT}, dt::FT) where FT
    nz    = grid.nz
    dz    = FT(grid.dz)
    Z_max = FT(nz) * dz

    l   = zeros(FT, nz)
    K_h = zeros(FT, nz)
    mixing_length_blackadar!(l, grid.centers_z, tke.l_inf)
    for k in 1:nz
        K_h[k] = tke.c_h * l[k] * sqrt(max(grid.states.e[k], tke.e_min))
    end

    for i in eachindex(droplets.X)
        droplets.cell_id[i] < 1 && continue   # already settled out

        z = droplets.z_loc[i]
        k = clamp(Int(floor(z / dz)) + 1, 1, nz)

        if k < nz
            z_kc  = (FT(k) - FT(0.5)) * dz  # centre of cell k
            frac  = clamp((z - z_kc) / dz, FT(0), FT(1))
            K_loc = (FT(1) - frac) * K_h[k] + frac * K_h[k + 1]
            dKdz  = (K_h[k + 1] - K_h[k]) / dz
        else
            K_loc = K_h[nz]
            dKdz  = FT(0)
        end

        z_new = z + dKdz * dt + sqrt(FT(2) * K_loc * dt) * randn(FT)

        if z_new < FT(0)
            z_new = -z_new
        elseif z_new > Z_max
            z_new = FT(2) * Z_max - z_new
        end

        droplets.z_loc[i]   = z_new
        droplets.cell_id[i] = Int(floor(z_new / dz)) + 1
    end
end
