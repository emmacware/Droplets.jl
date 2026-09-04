include("nonlocal.jl")
export my25_stability_functions, diffuse_fields!
export turb_timestep!, tke_update!, bott_mixing_length!
export deardorff_mixing_length!
export edmfx_mixing_length!
export mynn_stability_functions
export mixing_length!




"""
    turb_timestep!(::DynON, grid, tke, constants, dt, dynamics, turbdata)

Advance the TKE-closure turbulence scheme by one timestep. Computes the
mixing length (`tke.mixing_length_scheme`) and Mellor-Yamada stability
functions/eddy diffusivities (`K_h`, `K_m`, `K_e`) from the current TKE
field, runs turbulent droplet diffusion, implicitly diffuses the prognostic
fields (thermodynamic pair, wind, TKE) via [`diffuse_fields!`](@ref), records
the resulting vertical-transport TKE tendency, then updates TKE via
production/dissipation with [`tke_update!`](@ref).

"""
function turb_timestep!(::DynON,grid::scm_eulerian_arrays{FT}, tke::tke_settings{FT}, constants::Constants{FT}, dt::FT, dynamics::dynamic_settings,
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
    turbdata.shear_production .= 0
    turbdata.buoyancy_production .= 0

    e_prev = grid.states.e .+ 0

    T_from_theta!(grid.states.T_tmp,grid.states.θ, grid.states.P, grid.states.qv,constants)

    mixing_length!(tke.mixing_length_scheme, l, grid, tke, constants)

    # calculate_Kh_Km!(turbscheme, K_h, K_m,K_e, l, N2,S2,SM,SH,GN,GH,grid, tke, constants)

    for k in 1:nz
        my25_stability_functions(l,K_h,K_m,K_e,GH,GM,SM,SH, k,tke,grid,constants)
    end
    run_droplet_diffusion!(dynamics.droplet_diffusion_scheme, dynamics.turbulent_droplet_diffusion_on,
        l, droplets, grid, tke, dt, K_h, e_prev)

    diffuse_fields!(grid, tke, K_h,K_m,K_e,constants, dt,turbdata)

    turbdata.transport .= (grid.states.e .- e_prev) ./ dt

    tke_update!(l,SM,SH,GM,GH,grid, tke,dt, constants, turbdata)

end

"""
    turb_timestep!(::DynOFF, grid, tke, constants, dt, dynamics, turbdata)

Turbulence scheme disabled.
"""
function turb_timestep!(::DynOFF,grid::scm_eulerian_arrays{FT}, tke::tke_settings{FT}, constants::Constants{FT}, dt::FT, dynamics::dynamic_settings,
    turbdata::turbulence_data) where FT
end






"""
    diffuse_thermo!(::ThetalQtVar, grid, tke, K_h, dt, dz, nz, turbdata, theta_surf_flux, constants, ρ_arg)

Implicitly diffuse θl and qt (conserved under phase change), then
reconstruct θ from the diffused θl and the pre-diffusion `ql_tmp`, and qv as
`qt_tmp - ql_tmp`.
"""
function diffuse_thermo!(::ThetalQtVar, grid, tke, K_h, dt, dz, nz, turbdata, theta_surf_flux, constants, ρ_arg)
    implicit_diffuse!(grid.states.θl_tmp, K_h, dt, dz, nz,turbdata, sfc_flux = theta_surf_flux, ρ = ρ_arg)
    grid.states.θ .= θ_from_θl.(grid.states.P,grid.states.θl_tmp,grid.states.ql_tmp,grid.states.qv,constants)

    implicit_diffuse!(grid.states.qt_tmp, K_h, dt, dz, nz,turbdata, sfc_flux = tke.LHF / (constants.L*grid.states.ρ[1]), ρ = ρ_arg)
    grid.states.qv .= grid.states.qt_tmp .- grid.states.ql_tmp
end

"""
    diffuse_thermo!(::ThetaQvVar, grid, tke, K_h, dt, dz, nz, turbdata, theta_surf_flux, constants, ρ_arg)

Implicitly diffuse θ and qv directly, then recompute `qt_tmp` and `θl_tmp`
from the diffused θ/qv and the (unmodified) `ql_tmp`.
"""
function diffuse_thermo!(::ThetaQvVar, grid, tke, K_h, dt, dz, nz, turbdata, theta_surf_flux, constants, ρ_arg)
    implicit_diffuse!(grid.states.θ, K_h, dt, dz, nz,turbdata, sfc_flux = theta_surf_flux, ρ = ρ_arg)

    # grid.states.qv .+= grid.states.ql_tmp
    implicit_diffuse!(grid.states.qv, K_h, dt, dz, nz,turbdata, sfc_flux = tke.LHF / (constants.L*grid.states.ρ[1]), ρ = ρ_arg)
    # grid.states.qv .-= grid.states.ql_tmp
    grid.states.qt_tmp .= grid.states.qv .+ grid.states.ql_tmp
    grid.states.θl_tmp .= θl.(grid.states.P,grid.states.θ,grid.states.qt_tmp,grid.states.qv,constants)
end

"""
    diffuse_fields!(grid, tke, K_h, K_m, K_e, constants, dt, turbdata)

One implicit-diffusion timestep for the full prognostic state:
- `ThetalQtVar` or `ThetaQvVar`, with surface heat flux `tke.SHF` and
  latent heat flux `tke.LHF`
- horizontal wind components `grid.wind.u`/`grid.wind.v`, with surface
  momentum flux from `tke.u_star` and an explicit Coriolis update afterward
  using `tke.geostrophic_u`/`tke.geostrophic_v`
- TKE `grid.states.e`, with surface flux `2.5 * tke.u_star^3`, floored at 0

Diffusion is density-weighted (`ρ_arg = grid.states.ρ`) when
`tke.density_weighted_diffusion` is true, else unweighted. Density
`grid.states.ρ` is re-diagnosed after the thermodynamic update via
`ρ_calc_θ!`.
"""
function diffuse_fields!(grid,tke, K_h,K_m,K_e, constants, dt,turbdata)
    nz = grid.nz
    dz = grid.dz
    FT = eltype(grid.states.θ)
    # ρ = grid.states.ρ .+ 0
    spd = max(sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2), FT(1e-10))


    surface_zonal_momentum_flux = -grid.wind.u[1]*tke.u_star^2 /spd
    surface_meridional_momentum_flux = - grid.wind.v[1]*tke.u_star^2 / spd
    theta_surf_flux = tke.SHF / (grid.states.ρ[1] * constants.Cp_air * (grid.states.P[1]/constants.P0)^(constants.Rd / constants.Cp_air))

    ρ_arg = tke.density_weighted_diffusion ? grid.states.ρ : nothing

    diffuse_thermo!(tke.thermo_variable, grid, tke, K_h, dt, dz, nz, turbdata, theta_surf_flux, constants, ρ_arg)

    implicit_diffuse!(grid.wind.u,  K_m, dt, dz, nz,turbdata, sfc_flux = surface_zonal_momentum_flux, ρ = ρ_arg)
    implicit_diffuse!(grid.wind.v,  K_m, dt, dz, nz,turbdata, sfc_flux = surface_meridional_momentum_flux, ρ = ρ_arg)

    ρ_calc_θ!(grid.states.ρ,grid.states.P, grid.states.θ, grid.states.qv, constants)

    coriolis_parameter = 2 * constants.Ω * sind(31.5) #0.76e-4
    for k in 1:nz
        z = grid.centers_z[k]
        turbdata.du[k] = coriolis_parameter * (grid.wind.v[k] - tke.geostrophic_v(z))
        turbdata.dv[k] = -coriolis_parameter * (grid.wind.u[k] - tke.geostrophic_u(z))
    end

    grid.wind.u .+= turbdata.du .* dt
    grid.wind.v .+= turbdata.dv .* dt


    implicit_diffuse!(grid.states.e, K_e, dt, dz, nz,turbdata, sfc_flux=2.5 * tke.u_star^3, ρ = ρ_arg)
    grid.states.e .= max.(grid.states.e, 0.0)

end


"""
    deardorff_mixing_length!(l, grid, tke, constants)

Deardorff (1980)-style mixing length. 
"""
function deardorff_mixing_length!(l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    nz  = grid.nz
    dry = tke.dry_buoyancy
    c_n = 0.1
    for k in 1:nz
        N2 = calculate_buoyancy_frequency(grid, k, constants; dry)
        if N2 > 0
            l[k] = min(70, c_n * sqrt(grid.states.e[k] / N2))
        else
            l[k] = 70
        end
        lf = max(l[k], 1.0)
        l[k] = lf
    end
end

"""
    bott_mixing_length!(l, grid, tke, constants)

Bott (1996) mixing length
"""
function bott_mixing_length!(l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    nz  = grid.nz
    dz = grid.dz
    raw_inv = findfirst(k -> grid.states.qv[k] < 0.008, 1:nz)
    z_inv_idx = isnothing(raw_inv) ? nz : max(raw_inv - 1, 1)

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
        l[k] = k == 1 ? lf : max(lf, dz)
    end
end

"""
    edmfx_mixing_length!(l, grid, tke, constants)

adapted from ClimaAtmos's EDMF-X closure
(`ClimaAtmos.jl/src/prognostic_equations/edmfx_closures.jl`, `mixing_length`;
coefficients from Lopez-Gomez et al. 2020, https://doi.org/10.1029/2020MS002162).

Simplifications relative to the ClimaAtmos source (single-column, no EDMF updraft/
environment decomposition here):
  - l_W: simplified to neutral von Kármán scaling (κz), Blackadar (1962)-blended
    with a fixed asymptotic scale `tke.l_inf` so it doesn't grow unbounded with
    height; the Businger/Gryanik Monin-Obukhov stability correction was not ported.
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

        l_W = vk * z * tke.l_inf / (vk * z + tke.l_inf)

        N2 = calculate_buoyancy_frequency(grid, k, constants, dry=false)
        S2 = S2_flow_deformation(grid, k, constants, tke)
        sqrt_e = sqrt(e)
        a_pd = tke.c_m * (S2 - N2) * sqrt_e   # Pr_t = 1 simplification
        l_TKE = a_pd > eps(FT) ? sqrt(tke.c_d * e * sqrt_e / a_pd) : zero(FT)

        l_N = z
        if N2 > eps(FT) && e > eps(FT)
            l_N = min(sqrt(tke.c_b * e) / sqrt(N2), z)
        end

        l[k] = max(min(l_W, l_TKE > 0 ? l_TKE : l_W, l_N), min(dz, z))
    end
end

"""
    mixing_length!(scheme, l, grid, tke, constants)

Dispatch on `tke.mixing_length_scheme` to fill `l` in place:
`DeardorffMixing` → [`deardorff_mixing_length!`](@ref),
`BottMixing` → [`bott_mixing_length!`](@ref),
`EDMFXMixing` → [`edmfx_mixing_length!`](@ref).
"""
mixing_length!(::DeardorffMixing, l, grid, tke, constants) = deardorff_mixing_length!(l, grid, tke, constants)
mixing_length!(::BottMixing, l, grid, tke, constants) = bott_mixing_length!(l, grid, tke, constants)
mixing_length!(::EDMFXMixing, l, grid, tke, constants) = edmfx_mixing_length!(l, grid, tke, constants)




"""
    my25_stability_functions(l, K_h, K_m, K_e, GH, GM, SM, SH, k, tke, grid, constants)

Mellor-Yamada level-2.5 (1982) stability functions at level `k`, following
Bott (1996,1997).
"""
function my25_stability_functions(l::Vector{FT},K_h::Vector{FT}, K_m::Vector{FT},K_e::Vector{FT},GH,GM,SM,SH, k,tke,grid,constants) where FT
    compute_ql_at_cell!(grid.states, k,constants)
    # dry = tke.dry_buoyancy
    # Π = (grid.states.P[k] / constants.P0) ^ (constants.Rd / constants.Cp_air)

    θl_k = θl(grid.states.P[k],grid.states.T_tmp[k],grid.states.ql_tmp[k],grid.states.qv[k],constants)

    q      = sqrt(2 * max(grid.states.e[k], tke.e_min))
    # N2  = calculate_buoyancy_frequency(grid, k, constants, dry=false)
    S2  = S2_flow_deformation(grid, k, constants, tke)
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



"""
    tke_update!(l, SM, SH, GM, GH, grid, tke, dt, constants, turbdata)

Advance `grid.states.e` (TKE) one timestep at each level via a production/
dissipation update, implicit in the dissipation term
"""
function tke_update!(l,SM,SH,GM,GH,grid, tke,dt, constants, turbdata::turbulence_data)
    nz = grid.nz
    dz = grid.dz

    for k in 1:nz
        q         = sqrt(2 * grid.states.e[k])
        turbdata.shear_production[k]    = SM[k] * GM[k] * q^3 / l[k]
        turbdata.buoyancy_production[k] = SH[k] * GH[k] * q^3 / l[k]
        prod      = turbdata.shear_production[k] + turbdata.buoyancy_production[k]
        diss_rate = tke.my_diss * 2 * q / l[k]   # s⁻¹ = ε/e
        # grid.states.e[k] = max((grid.states.e[k] + dt * prod) / (1 + dt * diss_rate), tke.e_min)
        e_new = max((grid.states.e[k] + dt * prod) / (1 + dt * diss_rate), zero(eltype(l)))
        grid.states.eps[k] = diss_rate * e_new   # ε = (ε/e)*e [m²/s³]
        grid.states.e[k] = e_new
    end

    return nothing
end



"""
    bott1997term(grid, k, constants)

Buoyancy-flux term at level `k` following Bott (1997), in place of
a directly computed dry/moist buoyancy frequency.
"""
function bott1997term(grid, k, constants)
    nz  = grid.nz
    dz  = grid.dz

    θl_arr = grid.states.θl_tmp
    qv_arr = grid.states.qv
    ql_arr = grid.states.ql_tmp
    qt(i)  = qv_arr[i] + ql_arr[i]

    sm3(arr, i) = (arr(max(i-1,1)) + arr(i) + arr(min(i+1,nz))) / 3

    θl_k  = θl_arr[k]

    if k == 1
        # Second-order one-sided (3-point, forward) difference, actually centered at z_1
        dθldz_k = (-3*sm3(i->θl_arr[i],1) + 4*sm3(i->θl_arr[i],2) - sm3(i->θl_arr[i],3)) / (2*dz)
        dqtdz_k = (-3*sm3(qt,1) + 4*sm3(qt,2) - sm3(qt,3)) / (2*dz)
    elseif k == nz
        dθldz_k = (3*sm3(i->θl_arr[i],nz) - 4*sm3(i->θl_arr[i],nz-1) + sm3(i->θl_arr[i],nz-2)) / (2*dz)
        dqtdz_k = (3*sm3(qt,nz) - 4*sm3(qt,nz-1) + sm3(qt,nz-2)) / (2*dz)
    else
        dθldz_k = (sm3(i->θl_arr[i],k+1) - sm3(i->θl_arr[i],k-1)) / (2*dz)
        dqtdz_k = (sm3(qt,k+1) - sm3(qt,k-1)) / (2*dz)
    end

    qt_k = qt(k)
    ql   = grid.states.ql_tmp[k]
    T    = grid.states.T_tmp[k]
    θ    = grid.states.θ[k]
    qsat = saturation_specific_humidity(θ, grid.states.P[k],grid.states.qv[k], constants)
    S_env = sat(grid.states.qv[k], grid.states.P[k]) / esat(T)

    # α  = ql > FT(1e-5) ? exp(0.6*(min(S_env*100,100)-100)) : FT(0)
    α  = exp(0.6*(min(S_env*100,100)-100))
    δ  = virtual_temp_coeff(constants)  # Rv/Rd - 1 (≈0.61), shared with T_virtual/θ_virtual in conversions.jl
    b1 = 1 + δ*qt_k - (δ+1)*ql
    a1 = (1 + constants.L^2 * qsat / (constants.Cp_air * constants.Rv * T^2))^(-1)
    a2 = a1 * constants.L * qsat / (constants.Rv * T * θ)
    b2 = (1 + δ*qt_k - 2*(δ+1)*ql) * (constants.L * θ) / (constants.Cp_air * T) - (δ+1)*θl_k
    b3 = δ * (θl_k + constants.L * θ * ql / (constants.Cp_air * T))

    return (b1 - b2*a2*α)*dθldz_k + (b3 + b2*a1*α)*dqtdz_k
end

"""
    saturation_specific_humidity(θ, P, q_vap, constants)

Saturation specific humidity at temperature and pressure P

"""
function saturation_specific_humidity(θ, P,q_vap, constants)
    T = T_from_theta(θ, P,q_vap, constants)
    es = esat(T) # saturation vapor pressure in Pa
    qsat = constants.ϵ*es/(P - (1-constants.ϵ)*es)
    return qsat
end

