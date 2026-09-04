export S2_flow_deformation
export calculate_buoyancy_frequency
export implicit_diffuse!
export tke_settings
export turbulent_droplet_diffusion!, turbulent_droplet_diffusion_wellmixed!, turbulent_droplet_diffusion_visser!
export dθldz
export run_droplet_diffusion!


"""
    tke_settings{FT}(; kwargs...)

Configuration and tunable coefficients for the prognostic-TKE turbulence
closure: surface forcing (`u_star`, `SHF`, `LHF`), geostrophic wind profiles
(`geostrophic_u`, `geostrophic_v`), Mellor-Yamada/Galperin coefficients
(`c_n`, `vk`, `my_diss`, `GH_lims`), mixing-length scheme coefficients for
`bott_mixing_length!` (`bott_α`, `bott_β`) and `edmfx_mixing_length!`
(`c_m`, `c_d`, `c_b`, `l_inf`), which mixing-length scheme and thermodynamic
variable pair to dispatch on (`mixing_length_scheme`, `thermo_variable`), and
knobs for the droplet-diffusion schemes in this file (`average_e_l_3pt`,
`droplet_diffusion_length_dz`, `tau_max`) and for
[`implicit_diffuse!`](@ref) (`density_weighted_diffusion`).
"""
struct tke_settings{FT<:AbstractFloat, GU, GV, Mix<:MixingLengthScheme, TVar<:ThermoVariable}
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
    l_inf::FT  # Blackadar (1962) asymptotic mixing length, caps edmfx_mixing_length!'s l_W = vk*z term
    mixing_length_scheme::Mix  # which of deardorff_/bott_/edmfx_mixing_length! computes l
    thermo_variable::TVar  # whether diffuse_fields! diffuses (θ,qv) directly or (θ_l,qt)
    average_e_l_3pt::Bool  # droplet diffusion: average e (and l, when used) over [z-1,z,z+1] vs. cell z alone
    droplet_diffusion_length_dz::Bool  # droplet diffusion: use dz as the length scale instead of the diagnosed mixing length l
    density_weighted_diffusion::Bool  # mass-weight implicit_diffuse! (ρK face diffusivity, dt/(ρdz²) per-cell) instead of plain ∂/∂z(K∂ϕ/∂z)
    tau_max::FT  # cap on turbulent_droplet_diffusion_wellmixed!'s OU decorrelation timescale
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
        GH_lims::Tuple{FT,FT}   = (FT(-0.4), FT(0.03)),
        c_m::FT                 = FT(0.14),
        c_d::FT                 = FT(0.22),
        c_b::FT                 = FT(0.4),
        l_inf::FT               = FT(70.0),
        mixing_length_scheme::MixingLengthScheme = EDMFXMixing(),
        thermo_variable::ThermoVariable          = ThetalQtVar(),
        average_e_l_3pt::Bool    = true,
        droplet_diffusion_length_dz::Bool = true,
        density_weighted_diffusion::Bool = false,
        tau_max::FT              = FT(120.0),
    ) where {FT<:AbstractFloat}
    tke_settings{FT, typeof(geostrophic_u), typeof(geostrophic_v), typeof(mixing_length_scheme), typeof(thermo_variable)}(
        e_min, u_star, SHF, LHF, geostrophic_u, geostrophic_v,
        dry_buoyancy, c_n, vk, bott_α, bott_β, my_diss, GH_lims, c_m, c_d, c_b, l_inf,
        mixing_length_scheme, thermo_variable, average_e_l_3pt, droplet_diffusion_length_dz,
        density_weighted_diffusion, tau_max)
end


"""
    S2_flow_deformation(grid, k, constants, tke)

Squared vertical wind shear |dU/dz|² at level `k`, used as the shear-production
term in the Mellor-Yamada stability functions. 
"""
function S2_flow_deformation(grid, k, constants, tke)
    if k == 1
        z1 = grid.centers_z[1]
        S2 = (tke.u_star / (tke.vk * z1))^2
        return S2
    elseif k == grid.nz
        du_dz = (3*grid.wind.u[end] - 4*grid.wind.u[end-1] + grid.wind.u[end-2]) / (2*grid.dz)
        dv_dz = (3*grid.wind.v[end] - 4*grid.wind.v[end-1] + grid.wind.v[end-2]) / (2*grid.dz)
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

"""
    calculate_buoyancy_frequency(grid, k, constants; dry=false)

Brunt-Väisälä (buoyancy) frequency squared, N², at level `k`. With
`dry=true`, uses a centered difference of dry potential temperature `θ`;
otherwise uses moist `θl_tmp` scaled by [`bott1997term`](@ref) (Bott
1997's liquid-water correction to the moist static stability).
"""
function calculate_buoyancy_frequency(grid, k, constants; dry::Bool=false)
    g  = constants.gconst

    if dry
        nz = grid.nz
        dz = grid.dz
        kp = min(k + 1, nz)
        km = max(k - 1, 1)
        dz_buo = (k == 1 || k == nz) ? dz : 2.0 * dz
        θ = grid.states.θ
        return (g / θ[k]) * (θ[kp] - θ[km]) / dz_buo
    end

    return (g / grid.states.θl_tmp[k]) * bott1997term(grid, k, constants)
end


"""
    dθldz(θl_t, grid, k, constants)

Centered (one-sided at the boundaries) vertical gradient of liquid-water
potential temperature `θl_t` at level `k`.
"""
function dθldz(θl_t, grid, k, constants)
    nz = grid.nz
    dz = grid.dz
    kp = min(k + 1, nz)
    km = max(k - 1, 1)
    dz_buo = (k == 1 || k == nz) ? dz : 2.0 * dz


    return (θl_t[kp] - θl_t[km]) / dz_buo

end





"""
    implicit_diffuse!(ϕ, K_centers, dt, dz, nz, turbdata; sfc_flux=nothing, ρ=nothing)

Advance the diffusion equation ∂ϕ/∂t = ∂/∂z(K ∂ϕ/∂z) (or the mass-weighted
∂/∂z(ρK ∂ϕ/∂z) form when `ρ` is given) by one implicit (Crank-Nicolson)
timestep in place, via the Thomas algorithm. Both boundaries are no-flux
unless `sfc_flux` is given.
"""
function implicit_diffuse!(ϕ::Vector{FT}, K_centers::Vector{FT},
    dt::FT, dz::FT, nz::Int,turbdata::turbulence_data;
    sfc_flux::Union{FT,Nothing} = nothing,
    ρ::Union{Vector{FT},Nothing} = nothing) where FT

    #zero out the arrays in turbulence_data
    turbdata.K_faces .= 0
    turbdata.rhs .= 0
    turbdata.a .= 0
    turbdata.b .= 0
    turbdata.c .= 0
    turbdata.d .= 0

    # Face diffusivities, mass-weighted (ρK, face-averaged) when ρ is given
    K_face = turbdata.K_faces

    if ρ === nothing
        for k in 2:nz
            K_face[k] = FT(0.5) * (K_centers[k-1] + K_centers[k])
        end
    else
        for k in 2:nz
            K_face[k] = FT(0.5) * (ρ[k-1]*K_centers[k-1] + ρ[k]*K_centers[k])
        end
    end
    K_face[1] = 0 # (bottom NoFlux)
    K_face[nz+1] = 0 # (top NoFlux)

    r_cell = ρ === nothing ? fill(dt / dz^2, nz) : dt ./ (ρ .* dz^2)

    rhs = turbdata.rhs#zeros(FT, nz)   # right-hand side (copy of ϕ)
    rhs .= ϕ

    a = turbdata.a#zeros(FT, nz)   # sub-diagonal
    b = turbdata.b#zeros(FT, nz)   # main diagonal
    c = turbdata.c#zeros(FT, nz)   # super-diagonal
    d = turbdata.d#zeros(FT, nz)   # RHS

    for k in 2:nz-1
        Km = K_face[k]
        Kp = K_face[k+1]
        r = r_cell[k]
        a[k] = -r/2 * Km
        b[k] = 1 + r/2 * (Km + Kp)
        c[k] = -r/2 * Kp
        d[k] = ϕ[k] + r/2 * (Kp*(ϕ[k+1]-ϕ[k]) - Km*(ϕ[k]-ϕ[k-1]))
    end

    r1 = r_cell[1]
    a[1] = 0               # no sub-diagonal for first cell
    b[1] = 1 + r1/2*(K_face[1]+K_face[2])
    c[1] = -r1/2 * K_face[2]
    d[1] = ϕ[1] + r1/2 * (K_face[2]*(ϕ[2]-ϕ[1]))
    rnz = r_cell[nz]
    a[nz] = -rnz/2 * K_face[nz]
    b[nz] = 1 + rnz/2 * K_face[nz]
    c[nz] = 0
    d[nz] = ϕ[nz] - rnz/2 * K_face[nz] * (ϕ[nz] - ϕ[nz-1])
    if sfc_flux !== nothing
        # sfc_flux is already a kinematic (density-normalized) flux at every call site
        # in this codebase (e.g. theta_surf_flux = SHF/(ρ[1]*Cp*Π), 2.5*u_star^3 for
        # TKE) -- dividing by ρ[1] here too would double-count the density
        # normalization, so this stays ρ-independent even when ρ is given for the
        # interior operator.
        d[1] += dt * (sfc_flux / dz)
    end


    c_star = zeros(FT, nz)
    d_star = zeros(FT, nz)
    c_star[1] = c[1] / b[1]
    d_star[1] = d[1] / b[1]
    bet = b[1]

    for k in 2:nz
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

"""
    turbulent_droplet_diffusion!(::DynOFF, l, droplets, grid, tke, dt)
    turbulent_droplet_diffusion!(::DynON, l, droplets, grid, tke, dt)

Stochastic sub-grid vertical transport of droplets via droplet specific `w_prime` (adapted from
Abade & Grabowski 2018, restricted to depend on TKE alone, with the
Gillespie (1996) two-normal weak-second-order integration scheme) Droplets at or above the
diagnosed inversion (`grid.states.qv < 0.008`, rediagnosed every time step) have `w_prime` clamped to
`≤ 0` to not coast in turbulence free zones
"""
function turbulent_droplet_diffusion!(::DynOFF,l,droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
end
function turbulent_droplet_diffusion!(::DynON,l,droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
    #edited from grabowski and abade 2018 to be just dependent on TKE
    #also using gillespie solution for dw'
    nz    = grid.nz
    dz    = grid.dz
    
    Ce = FT(0.63)
    Ct = FT(0.89)
    inv_idx = findfirst(k -> grid.states.qv[k] < 0.008, 1:grid.nz)
    inv_idx = isnothing(inv_idx) ? nz : inv_idx
    # compute_ql_at_cell!.(grid.states, 1:nz,constants)
    # grid.states.qv .+= grid.states.ql_tmp
    # inv_idx = findfirst(k -> grid.states.e[k] < tke.e_min, 1:grid.nz)
    z_inv   = isnothing(inv_idx) ? FT(Inf) : FT((inv_idx-1) * grid.dz)
    for z in 1:nz#inv_idx
        isempty(droplets.grid_range[z]) && continue
        k = droplets.I[droplets.grid_range[z]]

        if !isnothing(inv_idx) && z >= inv_idx
            droplets.w_prime[k] .= min.(droplets.w_prime[k], 0)
            continue
        end
        zm = max(z - 1, 1)
        zp = min(z + 1, nz)

        e_eff = tke.average_e_l_3pt ? (grid.states.e[z] + grid.states.e[zm] + grid.states.e[zp]) / FT(3) : grid.states.e[z]
        # N2 = calculate_buoyancy_frequency(grid, z, constants)
        term = 1
        # N2 > 0 && (term = 1 / (1 + 400 * N2))
        l_z = tke.droplet_diffusion_length_dz ? dz : (tke.average_e_l_3pt ? (l[z] + l[zm] + l[zp]) / FT(3) : l[z])
        tau   = Ct*l_z * sqrt(Ce / e_eff) * term
        sigma2 = (e_eff * FT(2/3))
        noise_amp = sqrt(2*sigma2 / tau)
        w_corr    = 1 - dt / (2*tau)
        dt_sqrt   = sqrt(dt)
        dt_32     = dt * dt_sqrt

        n1 = randn(FT, length(k))
        n2 = randn(FT, length(k))
        w0 = copy(droplets.w_prime[k])

        droplets.z_loc[k]   .+= w0 .* (dt * w_corr) .+ noise_amp .* FT(0.5) .* (n1 .+ n2 ./ sqrt(FT(3))) .* dt_32
        droplets.w_prime[k] .+= (-w0 .* (dt / tau) .+ noise_amp .* n1 .* dt_sqrt) .* w_corr
        # clamp!(droplets.w_prime[k], -2, 2)

        # for ki in k
        #     if droplets.z_loc[ki] > z_inv
        #             droplets.z_loc[ki] = 2*z_inv - droplets.z_loc[ki]
        #             droplets.w_prime[ki] = -droplets.w_prime[ki]
        #     elseif droplets.z_loc[ki] < 0
        #         droplets.z_loc[ki] = rand()* dz
        #         droplets.w_prime[ki] = -droplets.w_prime[ki]
        #     end
        # end

        droplets.cell_id[k] = clamp.(floor.(Int, droplets.z_loc[k] / grid.dz) .+ 1, -1, nz)
    end

    # compute_ql_at_cell!.(grid.states, 1:nz,constants)
    # grid.states.qv .= grid.states.qt_tmp .- grid.states.ql_tmp
    return
end

"""
    turbulent_droplet_diffusion_wellmixed!(::DynOFF, l, droplets, grid, tke, dt)
    turbulent_droplet_diffusion_wellmixed!(::DynON, l, droplets, grid, tke, dt)

Same Ornstein-Uhlenbeck droplet-transport parameterization as
[`turbulent_droplet_diffusion!`](@ref) (`tau = Ct*l*sqrt(Ce/e)`, capped at
`tke.tau_max`), but with an added well-mixed drift term
`(2/3) ∂e/∂z` (per Bahlali, Henry & Carissimo 2020) that
`turbulent_droplet_diffusion!` omits, so a uniform droplet distribution
stays uniform under a spatially varying `e`.
"""
function turbulent_droplet_diffusion_wellmixed!(::DynOFF,l,droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
end
function turbulent_droplet_diffusion_wellmixed!(::DynON,l,droplets::droplet_attributes_1d{FT},
    grid, tke::tke_settings{FT}, dt::FT) where FT
    nz = grid.nz
    dz = grid.dz
    Ce = FT(0.63)
    Ct = FT(0.89)
    Z_max = FT(nz) * dz

    for z in 1:nz
        isempty(droplets.grid_range[z]) && continue
        k = droplets.I[droplets.grid_range[z]]

        zm = max(z - 1, 1)
        zp = min(z + 1, nz)

        e_z   = tke.average_e_l_3pt ? (grid.states.e[z] + grid.states.e[zm] + grid.states.e[zp]) / FT(3) : grid.states.e[z]
        e_z < tke.e_min && (e_z = tke.e_min)
        l_z   = tke.droplet_diffusion_length_dz ? dz : (tke.average_e_l_3pt ? (l[z] + l[zm] + l[zp]) / FT(3) : l[z])

        tau    = min(Ct * l_z * sqrt(Ce / e_z), tke.tau_max)
        sigma2 = FT(2/3) * e_z
        b  = sqrt(2 * sigma2 / tau)
        w_corr  = 1 - dt / (2*tau)
        dt_sqrt = sqrt(dt)
        dt_32   = dt * dt_sqrt

        # ∂sigma2/∂z = (2/3) ∂e/∂z, centered on the same zm/zp used above
        dz_eff = (z == 1 || z == nz) ? dz : FT(2)*dz
        drift  = FT(2/3) * (grid.states.e[zp] - grid.states.e[zm]) / dz_eff * dt

        n1 = randn(FT, length(k))
        n2 = randn(FT, length(k))
        w0 = copy(droplets.w_prime[k])

        droplets.z_loc[k]   .+= w0 .* (dt * w_corr) .+ drift .* dt .+ b .* FT(0.5) .* (n1 .+ n2 ./ sqrt(FT(3))) .* dt_32
        droplets.w_prime[k] .+= (-w0 .* (dt / tau) .+ b .* n1 .* dt_sqrt) .* w_corr .+ drift

        for ki in k
            if droplets.z_loc[ki] > Z_max
                droplets.z_loc[ki] = 2*Z_max - droplets.z_loc[ki]
                droplets.w_prime[ki] = -droplets.w_prime[ki]
            elseif droplets.z_loc[ki] < 0
                droplets.z_loc[ki] = rand()* dz
                droplets.w_prime[ki] = -droplets.w_prime[ki]
            end
        end

        droplets.cell_id[k] = clamp.(floor.(Int, droplets.z_loc[k] / grid.dz) .+ 1, 0, nz)
    end

    return
end

"""
    turbulent_droplet_diffusion_visser!(::DynOFF, droplets, grid, K_centers, dt)
    turbulent_droplet_diffusion_visser!(::DynON, droplets, grid, K_centers, dt)

Random-walk droplet transport directly from an eddy diffusivity profile
`K_centers` (Visser 1997): each droplet is displaced by a deterministic
drift `dK/dz * dt` plus Gaussian noise `sqrt(2*K_half*dt)`, 
"""
function turbulent_droplet_diffusion_visser!(::DynOFF, droplets::droplet_attributes_1d{FT},
    grid, K_centers::Vector{FT}, dt::FT) where FT
end
function turbulent_droplet_diffusion_visser!(::DynON, droplets::droplet_attributes_1d{FT},
    grid, K_centers::Vector{FT}, dt::FT) where FT
    # Visser, A.W. (1997), "Using random walk models to simulate the vertical
    # distribution of particles in a turbulent water column", Mar. Ecol. Prog. Ser.
    # 158, 275-281.
    #
    nz = grid.nz
    dz = grid.dz
    zc = grid.centers_z
    Z_max = FT(nz) * dz

    for z in 1:nz
        isempty(droplets.grid_range[z]) && continue
        k = droplets.I[droplets.grid_range[z]]

        zm = max(z - 1, 1)
        zp = min(z + 1, nz)
        dz_eff = (z == 1 || z == nz) ? dz : FT(2)*dz
        dKdz   = (K_centers[zp] - K_centers[zm]) / dz_eff

        # K evaluated at the half-step position, linearly interpolated between
        # cell centers (K_centers is only known at cell centers, not continuously)
        z_half = zc[z] + FT(0.5) * dKdz * dt
        K_half = if z_half <= zc[1]
            K_centers[1]
        elseif z_half >= zc[nz]
            K_centers[nz]
        else
            jlo  = clamp(floor(Int, (z_half - zc[1]) / dz) + 1, 1, nz - 1)
            frac = (z_half - zc[jlo]) / dz
            K_centers[jlo] + frac * (K_centers[jlo+1] - K_centers[jlo])
        end
        K_half = max(K_half, zero(FT))

        R  = randn(FT, length(k))
        Δz = dKdz * dt .+ R .* sqrt(2 * K_half * dt)

        droplets.z_loc[k] .+= Δz

        for ki in k
            if droplets.z_loc[ki] > Z_max
                droplets.z_loc[ki] = 2*Z_max - droplets.z_loc[ki]
            elseif droplets.z_loc[ki] < 0
                droplets.z_loc[ki] = -droplets.z_loc[ki]
            end
        end

        droplets.cell_id[k] = clamp.(floor.(Int, droplets.z_loc[k] / dz) .+ 1, 1, nz)
    end

    return
end

"""
    run_droplet_diffusion!(scheme, onoff, l, droplets, grid, tke, dt, K_h, e_prev)

Dispatch on `dynamics.droplet_diffusion_scheme` to the corresponding
turbulent droplet transport implementation in this file:
`NoDropletDiffusion` → no-op,
`OUDropletDiffusion` → [`turbulent_droplet_diffusion!`](@ref),
`WellMixedDropletDiffusion` → [`turbulent_droplet_diffusion_wellmixed!`](@ref),
`WeilDropletDiffusion` → `weil_turbulent_droplet_diffusion!`,
`VisserDropletDiffusion` → [`turbulent_droplet_diffusion_visser!`](@ref).
`onoff` is `dynamics.turbulent_droplet_diffusion_on`, an independent
`DynON`/`DynOFF` master switch forwarded to whichever scheme is picked.
"""
run_droplet_diffusion!(::NoDropletDiffusion, onoff, l, droplets, grid, tke, dt, K_h, e_prev) = nothing
run_droplet_diffusion!(::OUDropletDiffusion, onoff, l, droplets, grid, tke, dt, K_h, e_prev) =
    turbulent_droplet_diffusion!(onoff, l, droplets, grid, tke, dt)
run_droplet_diffusion!(::WellMixedDropletDiffusion, onoff, l, droplets, grid, tke, dt, K_h, e_prev) =
    turbulent_droplet_diffusion_wellmixed!(onoff, l, droplets, grid, tke, dt)
run_droplet_diffusion!(::WeilDropletDiffusion, onoff, l, droplets, grid, tke, dt, K_h, e_prev) =
    weil_turbulent_droplet_diffusion!(onoff, l, droplets, grid, tke, dt, e_prev)
run_droplet_diffusion!(::VisserDropletDiffusion, onoff, l, droplets, grid, tke, dt, K_h, e_prev) =
    turbulent_droplet_diffusion_visser!(onoff, droplets, grid, K_h, dt)

"""
    stochastic_jump_diffusion!(::DynON, grid, droplets, K_centers, dt, dz, nz)

Discrete 1D Monte Carlo analogue of Eulerian K-diffusion: only cloud
droplets (radius < 40 μm) individually jump to the cell above/below with
probability proportional to the face diffusive flux `K_face*dt/dz²`
(renormalized if the two probabilities would sum above 1), independent of
each droplet's multiplicity `ξ`.
"""
function stochastic_jump_diffusion!(::DynON,grid,droplets::droplet_attributes_1d{FT}, K_centers::Vector{FT}, dt::FT, dz::FT, nz::Int) where FT
    # discrete 1D Monte Carlo analogue of Eulerian K-diffusion

    K_face = zeros(FT, nz + 1)
    for k in 2:nz
        K_face[k] = FT(0.5) * (K_centers[k-1] + K_centers[k])
    end
    K_face[1]    = 0  # no-flux bottom, matching implicit_diffuse!
    K_face[nz+1] = 0  # no-flux top, matching implicit_diffuse!

    for z in 1:nz
        isempty(droplets.grid_range[z]) && continue
        k_all = droplets.I[droplets.grid_range[z]]

        r = volume_to_radius.(droplets.X[k_all])
        cloud_mask = r .< FT(40e-6)   # only cloud droplets < 40 μm, matching tdroptests!
        k = k_all[cloud_mask]
        isempty(k) && continue

        # probability of crossing the upper/lower face this step, driven by the flux there
        p_up   = K_face[z+1] * dt / dz^2
        p_down = K_face[z]   * dt / dz^2
        p_tot  = p_up + p_down
        if p_tot > 1

            p_up   /= p_tot
            p_down /= p_tot
        end

        for ki in k
            u = rand(FT)
            if u < p_down && z > 1
                droplets.cell_id[ki] = z - 1
                droplets.z_loc[ki]   -=dz# (z - FT(1.5)) * dz
            elseif u < p_down + p_up && z < nz
                droplets.cell_id[ki] = z + 1
                droplets.z_loc[ki]   += dz# (z + FT(0.5)) * dz
            end
            # else: stays put, nothing to change
        end
    end

    return nothing
end

"""
    partmc_jump_diffusion!(::DynON, grid, droplets, K_centers, dt, dz, nz)

partmc-inspired stochastic diffusion of cloud droplets (radius < 40
μm) between adjacent cells. ... Differences due to weighting handling, 
not sure if this should be used yet
"""
function partmc_jump_diffusion!(::DynON,grid,droplets::droplet_attributes_1d{FT}, K_centers::Vector{FT}, dt::FT, dz::FT, nz::Int) where FT
    # partmc stochastic diffusion

    K_face = zeros(FT, nz + 1)
    for k in 2:nz
        K_face[k] = FT(0.5) * (K_centers[k-1] + K_centers[k])
    end
    K_face[1] = 0; K_face[nz+1] = 0  # no-flux top/bottom, matching implicit_diffuse!

    Ncell  = zeros(FT, nz)
    kcloud = Vector{Vector{Int}}(undef, nz)
    for z in 1:nz
        if isempty(droplets.grid_range[z])
            kcloud[z] = Int[]
            continue
        end
        k_all = droplets.I[droplets.grid_range[z]]
        r  = volume_to_radius.(droplets.X[k_all])
        kc = k_all[r .< FT(40e-6)]
        kcloud[z] = kc
        Ncell[z]  = isempty(kc) ? zero(FT) : FT(sum(droplets.ξ[kc]))
    end

    Fbar = zeros(FT, nz + 1)
    for z in 1:nz-1
        Fbar[z+1] = K_face[z+1] * dt / dz^2 * (Ncell[z] - Ncell[z+1])
    end

    for z in 1:nz
        isempty(kcloud[z]) && continue

        p_up   = (Fbar[z+1] > 0 && Ncell[z] > 0) ? min(Fbar[z+1] / Ncell[z], one(FT)) : zero(FT)
        p_down = (z > 1 && Fbar[z] < 0 && Ncell[z] > 0) ? min(-Fbar[z] / Ncell[z], one(FT)) : zero(FT)
        if p_up + p_down > 1
            p_up   /= (p_up + p_down)
            p_down /= (p_up + p_down)
        end

        for ki in kcloud[z]
            share = droplets.ξ[ki] / Ncell[z]
            u = rand(FT)
            if u < p_down * share && z > 1
                droplets.cell_id[ki] = z - 1
                droplets.z_loc[ki]   -=dz# (z - FT(1.5)) * dz
            elseif u < p_down * share + p_up * share && z < nz
                droplets.cell_id[ki] = z + 1
                droplets.z_loc[ki]   += dz#(z + FT(0.5)) * dz
            end
            # else: stays put, nothing to change
        end
    end

    return nothing
end
