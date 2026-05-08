include("nonlocal.jl")
export implicit_diffuse!, turbulent_droplet_diffusion!
export calculate_buoyancy_frequency, mixing_length!
export duynkerke_tke_timestep!, calculate_Ri, S2_flow_deformation
export calculate_Pr
export find_bl_height
export my25_stability_functions
export Pr_t_func,Ri, X
export AbstractTurbulencScheme, ProgTKE, MellorYamada, SmagorinskyLilly
export AbstractMixingLengthFunction, Deardorff, SL, Nonlocal
export turb_timestep!, Blackadar, KNMI, KNMIturb, tke_update!, apply_entrainment!






function turb_timestep!(turbscheme::AbstractTurbulenceScheme, grid::scm_eulerian_arrays{FT}, tke::tke_settings{FT}, constants::Constants{FT}, dt::FT, scmsettings) where FT
    nz = grid.nz
    dz = grid.dz
    droplets = grid.states.droplets
    K_m = zeros(FT, nz)
    K_h = zeros(FT, nz)
    l = zeros(FT, nz)
    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, 1:nz)
    mixing_length!(tke.mixing_length_function, l, grid, tke, constants)

    calculate_Kh_Km!(turbscheme, K_h, K_m, l, grid, tke, constants)

    if scmsettings.rho_weighted_diffusion
        diffuse_ρ_fields!(grid, tke, K_h, K_m, constants, dt)
    else
        diffuse_fields!(grid, tke, K_h, K_m, constants, dt)
    end

    # apply_counter_gradient!(grid, tke, K_h, constants, dt)
    tke_update!(turbscheme, K_m, K_h, l, grid, tke, dt, constants)

    if scmsettings.turbulent_droplet_diffusion_on
        turbulent_droplet_diffusion!(l, droplets, grid, tke, dt)
    end
end




function mixing_length!(::Deardorff, l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    c_n = FT(0.76)
    nz  = grid.nz
    dry = tke.dry_buoyancy
    for k in 1:nz
        N2 = calculate_buoyancy_frequency(grid, k, constants; dry)
        if N2 > 0
            l[k] = min(tke.l_inf, c_n * sqrt(grid.states.e[k] / N2))
        else
            l[k] = tke.l_inf
        end
        l[k] = max(l[k], 1e-3)
    end
end

function mixing_length!(::SL, l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    dry = tke.dry_buoyancy
    for k in 1:grid.nz
        N2   = calculate_buoyancy_frequency(grid, k, constants; dry)
        S2   = S2_flow_deformation(grid, k, constants)
        Pr_t = Pr_t_func(S2, N2)
        l[k] = N2 > 0 ? 0.2 * grid.dz * max(0, 1 - N2/(Pr_t * 2 * S2))^(1/4) : 0.2 * grid.dz
    end
end

function mixing_length!(::Nonlocal, l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    vonkarman = 0.4
    nz = grid.nz
    for k in 1:nz
        z = grid.centers_z[k]
        len = 30 + 270 * (1 - exp(-z/100))
        lc = 1/(vonkarman * z) + 1/len
        l[k] = 1/lc
    end
end

function mixing_length!(::Blackadar, l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    vonkarman = 0.4
    nz = grid.nz
    for k in 1:nz
        z = grid.centers_z[k]
        l[k] = tke.l_inf / (1 + tke.l_inf/(vonkarman * z))
    end
end

function mixing_length!(::KNMI, l::Vector{FT}, grid, tke::tke_settings{FT}, constants) where FT<:AbstractFloat
    vonkarman = 0.4
    nz  = grid.nz
    c_e = 0.1
    dry = tke.dry_buoyancy
    for k in 1:nz
        z  = grid.centers_z[k]
        N2 = calculate_buoyancy_frequency(grid, k, constants; dry)
        one_over_l = 1/(vonkarman * z) + 1/300
        if N2 < 0
            l[k] = 1/one_over_l
        else
            one_over_l_new = 1/one_over_l + 1/(c_e*sqrt(grid.states.e[k]/N2))
            l[k] = 1/one_over_l_new
        end
    end
end


function calculate_Kh_Km!(::ProgTKE, K_h::Vector{FT}, K_m::Vector{FT}, l::Vector{FT}, grid, tke, constants) where FT<:AbstractFloat
    nz    = grid.nz
    delta = grid.dz
    e     = grid.states.e
    for k in 1:nz
        sqrte  = sqrt(max(e[k], tke.e_min))
        K_m[k] = tke.c_m * l[k] * sqrte
        K_h[k] = K_m[k] * (1 + 2 * l[k] / delta)  # Deardorff 1980
    end
end

function tke_update!(::ProgTKE, K_m,K_h,l,grid, tke,dt, constants)
    nz = grid.nz
    FT = eltype(grid.states.θ)
    dz = FT(grid.dz)

    e = grid.states.e
    ρ = grid.states.ρ .+ 0

    ρe = ρ .* e
    implicit_diffuse!(ρe, K_m, dt, dz, nz, sfc_flux = 2.5 * ρ[1] * tke.u_star^3)
    e .= ρe ./ grid.states.ρ
    e .= max.(e, tke.e_min)
    grid.states.e .= e

    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, 1:nz)

    dry = tke.dry_buoyancy
    for k in 1:grid.nz
        c_eps = 0.19 + 0.51 * (l[k] / grid.dz)
        S2   = S2_flow_deformation(grid, k, constants)
        N2   = calculate_buoyancy_frequency(grid, k, constants; dry)
        Ps   = K_m[k] * S2
        Pb   = -K_h[k] * N2
        diss = c_eps * sqrt(max(grid.states.e[k], 1e-9)) / l[k]
        grid.states.e[k] = (grid.states.e[k] + dt*(Ps + Pb)) / (1 + dt*diss)
        grid.states.e[k] = max(grid.states.e[k], 1e-9)
    end

    return nothing
end

function calculate_Kh_Km!(::KNMIturb, K_h::Vector{FT}, K_m::Vector{FT},l::Vector{FT}, grid, tke, constants) where FT<:AbstractFloat
    nz = grid.nz
    dz = grid.dz
    e  = grid.states.e
    c = 0.516

    Pr = FT(0.77)
    for k in 1:nz
        sqrte  = sqrt(e[k])
        K_m[k] = c * l[k] * sqrte
        K_h[k] = K_m[k] / Pr
    end
    # MOST neutral-limit floor at the surface cell: K ≥ κ·u*·z
    K_m[1] = max(K_m[1], FT(0.4) * tke.u_star * grid.centers_z[1])
    K_h[1] = max(K_h[1], K_m[1] / Pr)
end

function tke_update!(::KNMIturb, K_m,K_h,l,grid, tke,dt, constants)
    nz = grid.nz
    dz = grid.dz
    e = grid.states.e
    ρ = grid.states.ρ .+ 0


    implicit_diffuse!(grid.states.e, K_m, dt, dz, nz, sfc_flux = 2.5 * tke.u_star^3)
    grid.states.e .= max.(grid.states.e, tke.e_min)

    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, 1:nz)

    for k in 1:grid.nz
        S2   = S2_flow_deformation(grid, k, constants)
        Ps   = K_m[k] * S2
        Pb   = moist_buoyancy_production(grid, k, K_h[k], constants)
        diss = (grid.states.e[k])^(3/2)/(15*l[k])
        grid.states.e[k] = (grid.states.e[k] + dt*(Ps + Pb)) / (1 + dt*diss)
        grid.states.e[k] = max(grid.states.e[k], 1e-9)
    end

    return nothing
end



function diffuse_fields!(grid,tke, K_h, K_m, constants, dt)
    nz = grid.nz
    dz = grid.dz
    FT = eltype(grid.states.θ)
    # ρ = grid.states.ρ .+ 0

    surface_zonal_momentum_flux = -grid.wind.u[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    surface_meridional_momentum_flux = - grid.wind.v[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    theta_surf_flux = tke.SHF / (grid.states.ρ[1] * constants.Cp_air * (grid.states.P[1]/constants.P0)^(constants.Rd / constants.Cp_air))

    # Diffuse θl (conserved under phase changes) rather than dry θ.
    # In the cloud layer θ = θl + L*ql/(Cp*Π); diffusing dry θ moves latent heat out of
    # the cloud each step, erroneously eroding LWC. ql_tmp was updated in turb_timestep!.
    Π = (grid.states.P ./ constants.P0) .^ (constants.Rd / constants.Cp_air)
    θl = grid.states.θ .- constants.L .* grid.states.ql_tmp ./ (constants.Cp_air .* Π)
    implicit_diffuse!(θl, K_h, dt, dz, nz, sfc_flux = theta_surf_flux)
    grid.states.θ .= θl .+ constants.L .* grid.states.ql_tmp ./ (constants.Cp_air .* Π)

    implicit_diffuse!(grid.states.qv, K_h, dt, dz, nz, sfc_flux = tke.LHF / (constants.L*grid.states.ρ[1]))

    implicit_diffuse!(grid.wind.u,  K_m, dt, dz, nz, sfc_flux = surface_zonal_momentum_flux)
    implicit_diffuse!(grid.wind.v,  K_m, dt, dz, nz, sfc_flux = surface_meridional_momentum_flux)

    grid.states.ρ .= ρ_calc_θ(grid.states.P, grid.states.θ, grid.states.qv, constants)

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
end

function diffuse_ρ_fields!(grid,tke, K_h, K_m, constants, dt)
    nz = grid.nz
    dz = grid.dz
    FT = eltype(grid.states.θ)
    ρ = grid.states.ρ .+ 0

    surface_zonal_momentum_flux = -grid.wind.u[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)
    surface_meridional_momentum_flux = - grid.wind.v[1]*tke.u_star^2 / sqrt(grid.wind.u[1]^2 + grid.wind.v[1]^2)

    # # diffuse ρ-weighted quantities: ∂(ρϕ)/∂t = ∂/∂z(ρ K ∂ϕ/∂z)
    ρθ  = ρ .* grid.states.θ
    ρqv = ρ .* grid.states.qv
    ρu  = ρ .* grid.wind.u
    ρv  = ρ .* grid.wind.v
        # #using exner
    theta_surf_flux = tke.SHF / (constants.Cp_air * (grid.states.P[1]/constants.P0)^(constants.Rd / constants.Cp_air))


    ρqv_old = ρ .* grid.states.qv   # save before diffusion

    implicit_diffuse!(ρθ,  K_h, dt, dz, nz, sfc_flux = theta_surf_flux)
    implicit_diffuse!(ρqv, K_h, dt, dz, nz, sfc_flux = tke.LHF / constants.L)
    implicit_diffuse!(ρu,  K_m, dt, dz, nz, sfc_flux = ρ[1]*surface_zonal_momentum_flux)
    implicit_diffuse!(ρv,  K_m, dt, dz, nz, sfc_flux = ρ[1]*surface_meridional_momentum_flux)

    # update ρ from moisture change only
    ρ .+= ρqv .- ρqv_old

    # then divide back
    grid.states.θ  .= ρθ  ./ ρ
    grid.states.qv .= ρqv ./ ρ
    grid.wind.u    .= ρu  ./ ρ
    grid.wind.v    .= ρv  ./ ρ
    grid.states.ρ .= ρ_calc_θ(grid.states.P, grid.states.θ, grid.states.qv, constants)

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
end


# The formula implemented is from Li et al. (JAS 2015, DOI: 10.1175/JAS-D-14-0335.1, their Eq. 39),
# with a reformulation and correction of an algebraic error in their expression:
Ri(S2,N2) = N2 / max(2*S2, 1e-12)
X(S2,N2) = 0.98 + 4.076923076923077 * Ri(S2,N2) #Pr_n is 0.98 for now, ω_pr is 4.076923076923077 in ClimaParams
Pr_t_func(S2,N2) = (X(S2,N2) + sqrt(max(X(S2,N2)^2 - 4*0.98*Ri(S2,N2), 0))) / 2

function calculate_Kh_Km!(::SmagorinskyLilly, K_h::Vector{FT}, K_m::Vector{FT}, l_phys::Vector{FT}, grid, tke, constants) where FT<:AbstractFloat
    nz = grid.nz
    dz = FT(grid.dz)

    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, 1:nz)
    c_smag = tke.c_smag
    delta  = cbrt(tke.delta_h^2 * dz)   # geometric mean filter scale (Δx*Δy*Δz)^(1/3)

    dry = tke.dry_buoyancy
    for k in 1:nz
        N2     = calculate_buoyancy_frequency(grid, k, constants; dry)
        S2     = S2_flow_deformation(grid, k, constants)
        S_norm = sqrt(max(S2, FT(0)))
        Ri_val = Ri(S2, N2)
        Pr_t   = Pr_t_func(S_norm, N2)
        fb     = Ri_val ≤ 0 ? FT(1) : max(FT(0), 1 - N2 / Pr_t / 2 / max(S2, eps(FT)))^(FT(1/4))
        νt     = c_smag^2 * (delta * fb)^2 * S_norm
        K_m[k] = νt
        K_h[k] = νt / Pr_t
        l_phys[k] = c_smag * delta * fb
        grid.states.e[k] = max(FT(0.5) * νt * S_norm, tke.e_min)
    end
    return nothing
end

function tke_update!(::SmagorinskyLilly,K_m,K_h,l, grid, tke,dt, constants)
    # No prognostic TKE in Smagorinsky scheme, so nothing to update here
    return nothing
end





"""
    my25_stability_functions(GH::FT) -> (Sm, Sh)

Mellor-Yamada (1982) level-2 algebraic stability functions for use in MY2.5.
Reference: Mellor & Yamada (1982), Rev. Geophys. Space Phys. 20, 851–875.

GH = −(l/q)² · N²  (positive for unstable stratification, negative for stable).
Returns (Sm, Sh) where  Km = l·q·Sm,  Kh = l·q·Sh.
GH is clamped to [−0.28, 0.0233] following Galperin et al. (1988).
Neutral values: Sm ≈ 0.393, Sh ≈ 0.494.
"""
function my25_stability_functions(GH::FT) where FT
    A1 = FT(0.92); A2 = FT(0.74)
    B1 = FT(16.6);  C1 = FT(0.08)

    GH = clamp(GH, FT(-0.28), FT(0.0233))

    # MY82 level-2 algebraic result for Sh (depends only on GH)
    Sh = A2 * (1 - 6*A1/B1) / (1 - 3*A2*(1 + 6*A1/B1)*GH)

    # MY82 level-2 algebraic result for Sm
    Sm = (A1*(1 - 3*C1 - 6*A1/B1) - A1*A2*(18*A1 + 9*A2)*GH/B1) / (1 - 9*A1*A2*GH)

    return Sm, Sh
end


function calculate_Kh_Km!(::MellorYamada, K_h::Vector{FT}, K_m::Vector{FT},l::Vector{FT}, grid, tke, constants) where FT<:AbstractFloat

    droplets = grid.states.droplets
    B1 = FT(16.6)

    nz  = grid.nz
    dz  = FT(grid.dz)
    e   = grid.states.e
    ρ   = copy(grid.states.ρ)

    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, 1:nz)
    dry = tke.dry_buoyancy
    # MY2.5 eddy diffusivities from level-2 stability functions
    for k in 1:nz
        q      = sqrt(2 * max(e[k], tke.e_min))
        N2     = calculate_buoyancy_frequency(grid, k, constants; dry)
        GH     = -(l[k] / q)^2 * N2
        Sm, Sh = my25_stability_functions(GH)
        K_m[k] = l[k] * q * Sm
        K_h[k] = l[k] * q * Sh
    end

end

function tke_update!(::MellorYamada, K_m,K_h,l,grid, tke,dt, constants)
    nz = grid.nz
    dz = grid.dz
    e = grid.states.e
    ρ = grid.states.ρ .+ 0
    E1 = 1.8    # TKE transport coefficient: Kq = E1·l·q
    B1 = 16.6   # MY2.5 dissipation constant
    
    K_q = zeros(eltype(e), nz)
    for k in 1:nz
        q      = sqrt(2 * max(e[k], tke.e_min))
        K_q[k] = E1 * l[k] * q
    end

    # Transport ρe with TKE diffusivity Kq (MY82: Kq = E1·l·q ≠ Kh)
    # ρe = ρ .* e
    implicit_diffuse!(e, K_q, dt, dz, nz,sfc_flux=2.5 * ρ[1] * tke.u_star^3)
    grid.states.e .= max.(grid.states.e, tke.e_min)

    # Recompute ql after diffusion for correct N² in production terms
    grid.states.ql_tmp .= compute_ql_at_cell.(grid.states, 1:nz)

    # TKE source/sink: production explicit, dissipation semi-implicit
    # ε = q³/(B1·l),   ε/e = 2q/(B1·l)  →  e_new = (e + dt·(Ps+Pb)) / (1 + dt·2q/(B1·l))
    dry = tke.dry_buoyancy
    for k in 1:nz
        q    = sqrt(2 * max(grid.states.e[k], tke.e_min))
        N2   = calculate_buoyancy_frequency(grid, k, constants; dry)
        S2   = S2_flow_deformation(grid, k, constants)
        Ps   = K_m[k] * S2
        Pb   = -K_h[k] * N2
        # diss = 2 * q / (B1 * l[k])
        GH     = -(l[k] / q)^2 * N2
        GM = l[k]^2/q^2 * S2
        Sm, Sh = my25_stability_functions(GH)
        diss = (Sm * GM + Sh * GH - 1/16.6) * (2 * grid.states.e[k])^(3/2) / l[k]  
        grid.states.e[k] = (grid.states.e[k] + dt * (Ps + Pb)) / (1 + dt * diss)
    end

    grid.states.e .= max.(grid.states.e, tke.e_min)

    return nothing
end


# function turbulent_droplet_diffusion!(droplets::droplet_attributes_1d{FT},
#     grid, tke::tke_settings{FT}, dt::FT) where FT
#     #edited from grabowski and abade 2017 eqn 10, Dziekan 2019 uses this too.
#     #their version builds on w' from the previous timestep
#     nz    = grid.nz
#     l_z = grid.dz
#     Ct = FT(0.63)
#     for z in 1:nz
#         droplets.grid_range[z] === nothing && continue
#         k     = droplets.I[droplets.grid_range[z]]
#         e_z   = grid.states.e[z]
#         tau   = l_z / (2π)^(FT(1/3)) * sqrt(Ct / e_z)
#         @assert dt < 2.0 * tau
#         # tau = grid.dz / sqrt(e_z) # eddy turnover time based on grid scale and local TKE
#         sigma = sqrt(e_z * FT(2/3))
#         w_amp = sqrt(1 - exp(-2 * dt / tau)) * sigma
#         droplets.w_prime[k] .*= exp(-dt / tau)
#         droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))
#     end

#     return
# end


# function turbulent_droplet_diffusion!(l,droplets::droplet_attributes_1d{FT},
#     grid, tke::tke_settings{FT}, dt::FT) where FT
#     #edited from grabowski and abade 2017 eqn 10, Dziekan 2019 uses this too.
#     #their version builds on w' from the previous timestep
#     nz    = grid.nz
#     # l_z = grid.dz
#     Ct = FT(0.63)
#     C_e = FT(0.89)
#     for z in 1:nz
#         l_z = min(l[z], grid.dz)
#         droplets.grid_range[z] === nothing && continue
#         k     = droplets.I[droplets.grid_range[z]]
#         e_z   = grid.states.e[z]
#         tau   = Ct * (l_z * C_e/ e_z)^(2/3) *l_z
#         # @assert dt < 2.0 * tau
#         # tau = grid.dz / sqrt(e_z) # eddy turnover time based on grid scale and local TKE
#         sigma = sqrt(e_z * FT(2/3))
#         w_amp = sqrt(2 * sigma^2/tau)
#         droplets.w_prime[k] .*= -1/tau * dt
#         droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))* sqrt(dt)
#     end

#     return
# end


function turbulent_droplet_diffusion!(l, droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
    # Abade et al. 2018 eq. 26: OU process dw' = (-w'/τ + ½ dσ²/dz)dt + sqrt(2σ²/τ)dW
    # Thomson (1987) well-mixed correction: drift += ½ d(σ²)/dz · τ
    # Without it, droplets accumulate in low-TKE regions (spurious downward drift).
    nz  = grid.nz
    dz  = FT(grid.dz)
    C_τ = FT(0.63)
    c_e = FT(0.89)
    e   = grid.states.e
    for z in 1:nz
        l_z = l[z]
        droplets.grid_range[z] === nothing && continue
        k     = droplets.I[droplets.grid_range[z]]
        e_z   = max(e[z], tke.e_min)
        ε     = (e_z / c_e)^FT(1.5) / l_z
        tau   = C_τ * cbrt(l_z^2 / ε)
        sigma = sqrt(e_z * FT(2/3))
        decay = exp(-dt / tau)
        w_amp = sigma * sqrt(1 - decay^2)

        e_up = z < nz ? max(e[z+1], tke.e_min) : e_z
        e_dn = z >  1 ? max(e[z-1], tke.e_min) : e_z
        dz_eff  = (z == 1 || z == nz) ? dz : 2*dz
        dsigma2_dz = FT(2/3) * (e_up - e_dn) / dz_eff
        drift_correction = FT(0.5) * dsigma2_dz * tau

        droplets.w_prime[k] .*= decay
        droplets.w_prime[k] .+= drift_correction * dt
        droplets.w_prime[k] .+= w_amp .* randn(FT, length(k))
    end

    return
end
