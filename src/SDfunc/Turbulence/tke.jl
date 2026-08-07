include("nonlocal.jl")
export my25_stability_functions, diffuse_fields!
export turb_timestep!, tke_update!, bott_mixing_length!
export deardorff_mixing_length!
export edmfx_mixing_length!
export mynn_stability_functions




function turb_timestep!(::DynON,grid::scm_eulerian_arrays{FT}, tke::tke_settings{FT}, constants::Constants{FT}, dt::FT, scmsettings,
    turbdata::turbulence_data) where FT
    nz = grid.nz
    dz = grid.dz
    droplets = grid.states.droplets

    l = turbdata.l
    SM = turbdata.SM
    SH = turbdata.SH
    GM = turbdata.GM
    GH = turbdata.GH
    K_h = turbdata.K_h
    K_m = turbdata.K_m
    K_e = turbdata.K_e
    l .= 0
    SM .= 0
    SH .= 0
    GM .= 0
    GH .= 0
    K_h .= 0
    K_m .= 0
    K_e .= 0
    
    e_prev = grid.states.e .+ 0

    T_from_theta!(grid.states.T_tmp,grid.states.θ, grid.states.P, grid.states.qv,constants)
    compute_ql_at_cell!.(grid.states, 1:nz,constants)

    # deardorff_mixing_length!(l, grid, tke, constants)
    bott_mixing_length!(l, grid, tke, constants)
    # edmfx_mixing_length!(l, grid, tke, constants)

    # calculate_Kh_Km!(turbscheme, K_h, K_m,K_e, l, N2,S2,SM,SH,GN,GH,grid, tke, constants)

    for k in 1:nz
        my25_stability_functions(l,K_h,K_m,K_e,GH,GM,SM,SH, k,tke,grid,constants)
    end
    turbulent_droplet_diffusion!(scmsettings.turbulent_droplet_diffusion_on,l, droplets, grid, tke, dt)
    # stochastic_jump_diffusion!(scmsettings.turbulent_droplet_diffusion_on,grid,droplets, K_h, dt, dz, nz)
    # partmc_jump_diffusion!(scmsettings.turbulent_droplet_diffusion_on,grid,droplets, K_h, dt, dz, nz)


    diffuse_fields!(grid, tke, K_h,K_m,K_e,constants, dt,turbdata)

    # apply_counter_gradient!(grid, tke, K_h, constants, dt)
    tke_update!(l,SM,SH,GM,GH,grid, tke,dt, constants)


    # turbulent_droplet_diffusion!(scmsettings.turbulent_droplet_diffusion_on,l, droplets, grid, tke, dt)
    # weil_turbulent_droplet_diffusion!(scmsettings.turbulent_droplet_diffusion_on,l, droplets, grid, tke, dt,e_prev)
    # partmc_jump_diffusion!(scmsettings.turbulent_droplet_diffusion_on,grid,droplets, K_h, dt, dz, nz)
end

function turb_timestep!(::DynOFF,grid::scm_eulerian_arrays{FT}, tke::tke_settings{FT}, constants::Constants{FT}, dt::FT, scmsettings,
    turbdata::turbulence_data) where FT
end






function diffuse_fields!(grid,tke, K_h,K_m,K_e, constants, dt,turbdata)
    nz = grid.nz
    dz = grid.dz
    FT = eltype(grid.states.θ)
    # ρ = grid.states.ρ .+ 0
    spd = max(sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2), FT(1e-10))


    surface_zonal_momentum_flux = -grid.wind.u[1]*tke.u_star^2 /spd
    surface_meridional_momentum_flux = - grid.wind.v[1]*tke.u_star^2 / spd
    theta_surf_flux = tke.SHF / (grid.states.ρ[1] * constants.Cp_air * (grid.states.P[1]/constants.P0)^(constants.Rd / constants.Cp_air))

    T_from_theta!(grid.states.T_tmp, grid.states.θ, grid.states.P, grid.states.qv,constants) 
    grid.states.θ .= θl.(grid.states.P,grid.states.T_tmp,grid.states.ql_tmp,grid.states.qv,constants)
    implicit_diffuse!(grid.states.θ, K_h, dt, dz, nz,turbdata, sfc_flux = theta_surf_flux)
    grid.states.θ .= θ_from_θl.(grid.states.P,grid.states.θ,grid.states.ql_tmp,grid.states.qv,constants)
    T_from_theta!(grid.states.T_tmp, grid.states.θ, grid.states.P, grid.states.qv,constants) 



    grid.states.qv .+= grid.states.ql_tmp
    implicit_diffuse!(grid.states.qv, K_h, dt, dz, nz,turbdata, sfc_flux = tke.LHF / (constants.L*grid.states.ρ[1]))
    grid.states.qv .-= grid.states.ql_tmp
# implicit_diffuse!(θl_t, K_h, dt, dz, nz, sfc_flux = theta_surf_flux)
    # θ_my_bot!(θl_t,grid, constants) #

    implicit_diffuse!(grid.wind.u,  K_m, dt, dz, nz,turbdata, sfc_flux = surface_zonal_momentum_flux)
    implicit_diffuse!(grid.wind.v,  K_m, dt, dz, nz,turbdata, sfc_flux = surface_meridional_momentum_flux)

    ρ_calc_θ!(grid.states.ρ,grid.states.P, grid.states.θ, grid.states.qv, constants)

    coriolis_parameter = 2 * constants.Ω * sind(31.5) #0.76e-4
    for k in 1:nz
        z = grid.centers_z[k]
        turbdata.du[k] = coriolis_parameter * (grid.wind.v[k] - tke.geostrophic_v(z))
        turbdata.dv[k] = -coriolis_parameter * (grid.wind.u[k] - tke.geostrophic_u(z))
    end

    grid.wind.u .+= turbdata.du .* dt
    grid.wind.v .+= turbdata.dv .* dt


    implicit_diffuse!(grid.states.e, K_e, dt, dz, nz,turbdata, sfc_flux=2.5 * tke.u_star^3)
    grid.states.e .= max.(grid.states.e, 0.0)

end


function deardorff_mixing_length!(l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    nz  = grid.nz
    dry = tke.dry_buoyancy
    c_n = 0.1
    for k in 1:nz
        N2 = calculate_buoyancy_frequency(grid, k, constants; dry)
        if N2 > 0
            l[k] = min(70, c_n * sqrt(grid.states.e[k] / N2))
        else
            l[k] = 70#tke.l_inf
        end
        lf = max(l[k], 1.0)
        # lf1 = min(lf, grid.dz)
        l[k] = lf
    end
end

function bott_mixing_length!(l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    nz  = grid.nz
    dz = grid.dz
    z_inv_idx = findfirst(k -> grid.states.qv[k] < 0.008, 1:nz) - 1

    z_inv = grid.centers_z[z_inv_idx]
    # l0 = alpha * z_inv
    for k in 1:nz
        z = grid.centers_z[k]
        
        if z > z_inv
            l0 = dz
        else
            l0 = tke.bott_α*z_inv - (tke.bott_α*z_inv - dz) * exp((z-z_inv)/(tke.bott_β*dz))

        end
        lf = tke.vk * z * l0 / (tke.vk * z + l0)
        l[k] = max(lf, dz)
    end
end

"""
    edmfx_mixing_length!(l, grid, tke, constants)

Alternative to `bott_mixing_length!`, adapted from ClimaAtmos's EDMF-X closure
(`ClimaAtmos.jl/src/prognostic_equations/edmfx_closures.jl`, `mixing_length`;
coefficients from Lopez-Gomez et al. 2020, https://doi.org/10.1029/2020MS002162).

`bott_mixing_length!`'s asymptotic scale is `bott_α * z_inv` -- proportional to the
*currently diagnosed* inversion height. That creates a destabilizing feedback: if the
boundary layer starts to shrink, the mixing length shrinks with it, weakening the very
mixing that would otherwise resist further shrinking. None of the three scales blended
here reference a diagnosed BL-top height at all -- each is a purely local function of
the cell's own shear, buoyancy, and TKE -- so this closure can't develop that feedback
by construction.

Simplifications relative to the ClimaAtmos source (single-column, no EDMF updraft/
environment decomposition here):
  - l_W: simplified to neutral von Kármán scaling (κz); the Businger/Gryanik
    Monin-Obukhov stability correction was not ported.
  - l_TKE: the production-dissipation balance omits the EDMF exchange term
    (`ᶜtke_exch`), which has no counterpart outside EDMF's updraft/environment split.
  - Turbulent Prandtl number fixed at 1 in the buoyancy-production term.
  - Scales are combined with a hard minimum rather than ClimaAtmos's smooth minimum.
  - Floored (not capped) at `dz`, matching this codebase's existing convention.
"""
function edmfx_mixing_length!(l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    nz = grid.nz
    dz = grid.dz
    vk = tke.vk

    for k in 1:nz
        z = grid.centers_z[k]
        e = max(grid.states.e[k], zero(FT))

        # l_W: near-surface (neutral von Kármán) scale
        l_W = vk * z

        # l_TKE: local production/dissipation balance, a_pd*l² = c_d*e^(3/2)
        # (no tke_exch term -- that's an EDMF updraft/environment exchange concept
        # with no single-column counterpart)
        N2 = calculate_buoyancy_frequency(grid, k, constants, dry=false)
        S2 = S2_flow_deformation(grid, k, constants)
        sqrt_e = sqrt(e)
        a_pd = tke.c_m * (S2 - N2) * sqrt_e   # Pr_t = 1 simplification
        l_TKE = a_pd > eps(FT) ? sqrt(tke.c_d * e * sqrt_e / a_pd) : zero(FT)

        # l_N: buoyancy (static-stability) limited scale, capped by wall distance
        l_N = z
        if N2 > eps(FT) && e > eps(FT)
            l_N = min(sqrt(tke.c_b * e) / sqrt(N2), z)
        end

        l[k] = max(min(l_W, l_TKE > 0 ? l_TKE : l_W, l_N), dz)
    end
end




"""
    my25_stability_functions(GH::FT) -> (Sm, Sh)

Mellor-Yamada (1982)
Reference: Bott 96
"""
function my25_stability_functions(l::Vector{FT},K_h::Vector{FT}, K_m::Vector{FT},K_e::Vector{FT},GH,GM,SM,SH, k,tke,grid,constants) where FT
    compute_ql_at_cell!(grid.states, k,constants)
    # dry = tke.dry_buoyancy
    # Π = (grid.states.P[k] / constants.P0) ^ (constants.Rd / constants.Cp_air)

    θl_k = θl(grid.states.P[k],grid.states.T_tmp[k],grid.states.ql_tmp[k],grid.states.qv[k],constants)

    q      = sqrt(2 * max(grid.states.e[k], tke.e_min))
    N2  = calculate_buoyancy_frequency(grid, k, constants, dry=false)
    # N2 = moist_buoyancy_frequency(grid, k, constants)
    S2  = S2_flow_deformation(grid, k, constants)
    l2_2e = l[k]^2 / q^2

    # GH[k]  = - l2_2e * N2
    # GH[k]  = -constants.gconst * l2_2e * dθldz_k / θl
    GH[k] =  - constants.gconst * l2_2e / θl_k * bott1997term(grid,k,constants)
    GM[k]  = l2_2e * S2

    GH[k] = clamp(GH[k], tke.GH_lims[1], tke.GH_lims[2])
    GM[k] = min(GM[k],25.0*(0.03-GH[k]))
    A1 = FT(0.74); A2 = FT(-4.534); A3 = FT(0.902); A4 = FT(0.699); A5 = FT(-9.339)
    A6 = FT(-36.719); A7 = FT(187.441); A8 = FT(5.078);A9 = FT(-88.839)

    A10 = (1 + (A6 + A7*GH[k])*GH[k] + (A8 + A9*GH[k])*GM[k]) ^ (-1)

    SH[k] = (A1 + A2*GH[k] + A3*GM[k])*A10
    SM[k] = (A4 + A5*GH[k])*A10

    K_m[k] = l[k] * q * SM[k]
    K_h[k] = l[k] * q * SH[k]
    K_e[k] = 0.2 * l[k] * q

end


# function calculate_Kh_Km!(::MellorYamada, K_h::Vector{FT}, K_m::Vector{FT},K_e::Vector{FT},l::Vector{FT},N2,S2,SM,SH,GN,GH, grid, tke, constants) where FT<:AbstractFloat
#     for k in 1:nz
#         my25_stability_functions(GH,GM,SM,SH,k,tke,grid,constants)
#     end

# end

function tke_update!(l,SM,SH,GM,GH,grid, tke,dt, constants)
    nz = grid.nz
    dz = grid.dz


    #do this again?
    # for k in 1:nz
    #     my25_stability_functions()
    # end

    # CURRENT (buggy denominator): diss = ε [m²/s³], so dt*diss has units m²/s², not dimensionless.
    # for k in 1:nz
    #     de   = (SM[k] * GM[k] + SH[k] * GH[k]) * (2 * grid.states.e[k])^(3/2) / l[k]
    #     diss = tke.my_diss * (2 * grid.states.e[k])^(3/2) / l[k]
    #     grid.states.e[k] = (grid.states.e[k] + dt * (de)) / (1 + dt * diss)
    # end

    # grid.states.e .= min.(grid.states.e, 0.4)

    # FIXED (correct semi-implicit): diss_rate = ε/e = 2q/(B1*l) [s⁻¹], so dt*diss_rate is dimensionless.
    # Derivation: linearize ε^(n+1) ≈ (ε^n/e^n)*e^(n+1), solve → e = (e + dt*prod)/(1 + dt*diss_rate).
    # Reference: MY82 dissipation formula ε=q³/(B1*l); semi-implicit form in Durran (1999) §2.4.
    for k in 1:nz
        # q = sqrt(2 * max(grid.states.e[k], tke.e_min))
        q         = sqrt(2 * grid.states.e[k])
        prod      = (SM[k] * GM[k] + SH[k] * GH[k]) * q^3 / l[k]
        diss_rate = tke.my_diss * 2 * q / l[k]   # s⁻¹ = ε/e
        # grid.states.e[k] = max((grid.states.e[k] + dt * prod) / (1 + dt * diss_rate), tke.e_min)
        grid.states.e[k] = max((grid.states.e[k] + dt * prod) / (1 + dt * diss_rate), zero(eltype(l)))
    end
    grid.states.e .= max.(grid.states.e, tke.e_min)

    return nothing
end



function bott1997term(grid, k, constants)
    FT  = eltype(grid.states.θ)
    nz  = grid.nz
    dz  = grid.dz
    kp  = min(k + 1, nz)
    km  = max(k - 1, 1)
    dz_eff = (k == 1 || k == nz) ? dz : 2.0 * dz

    θl_k  = θl(grid.states.P[k],  grid.states.T_tmp[k],  grid.states.ql_tmp[k],grid.states.qv[k],  constants)
    θl_kp = θl(grid.states.P[kp], grid.states.T_tmp[kp], grid.states.ql_tmp[kp],grid.states.qv[kp], constants)
    θl_km = θl(grid.states.P[km], grid.states.T_tmp[km], grid.states.ql_tmp[km],grid.states.qv[km], constants)

    dθldz_k = (θl_kp - θl_km) / dz_eff
    qt_k  = grid.states.qv[k]  + grid.states.ql_tmp[k]
    qt_kp = grid.states.qv[kp] + grid.states.ql_tmp[kp]
    qt_km = grid.states.qv[km] + grid.states.ql_tmp[km]
    dqtdz_k = (qt_kp - qt_km) / dz_eff

    ql   = grid.states.ql_tmp[k]
    T    = grid.states.T_tmp[k]
    θ    = grid.states.θ[k]
    qsat = saturation_specific_humidity(θ, grid.states.P[k],grid.states.qv[k], constants)
    S_env = sat(grid.states.qv[k], grid.states.P[k]) / esat(T)

    # α  = ql > FT(1e-5) ? exp(0.6*(min(S_env*100,100)-100)) : FT(0)
    α  = exp(0.6*(min(S_env*100,100)-100))
    b1 = 1 + 0.61*qt_k - 1.61*ql
    a1 = (1 + constants.L^2 * qsat / (constants.Cp_air * constants.Rv * T^2))^(-1)
    a2 = a1 * constants.L * qsat / (constants.Rv * T * θ)
    b2 = (1 + 0.61*qt_k - 3.22*ql) * (constants.L * θ) / (constants.Cp_air * T) - 1.61*θl_k
    b3 = 0.61 * (θl_k + constants.L * θ * ql / (constants.Cp_air * T))

    return (b1 - b2*a2*α)*dθldz_k + (b3 + b2*a1*α)*dqtdz_k
end

function saturation_specific_humidity(θ, P,q_vap, constants)
    T = T_from_theta(θ, P,q_vap, constants)
    es = esat(T) # saturation vapor pressure in Pa
    # qsat = constants.ϵ * es / (P - es) # saturation specific humidity
    qsat = constants.ϵ*es/(P - (1-constants.ϵ)*es)
    return qsat
end

