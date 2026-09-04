export dP_dz, dρ_dry_dz_pysdm, init_droplets_scm

"""
dP_dz(P_z, ode_settings, z) computes the vertical pressure gradient dP/dz at a given height z, based on the current pressure P_z,
to solve for the initial hydrostatic pressure profile.
"""
function dP_dz(P_z, ode_settings, z)

    FT = eltype(P_z)
    constants,θl_func,qt_func,q_vap_density = ode_settings

    θ_z::FT = θl_func(z)
    q_vap::FT = qt_func(z)
    # q_vap_density = nothing (default): density uses q_vap
    # q_vap_density = something: use Single fixed reference q_vap for the density term across
    # the whole column instead, while still integrating the real θ(z)/q_vap(z) profile for
    # everything else
    q_vap_ρ::FT = q_vap_density === nothing ? q_vap : q_vap_density

    g::FT = constants.gconst

    T::FT = T_from_theta(θ_z,P_z,q_vap,constants)
    ρ::FT = ρ_ideal_gas(P_z, T, q_vap_ρ, constants)

    return -g * ρ
end


# recreating PySDM's kid1d density-profile construction
# (see PySDM/physics/hydrostatics/constant_g_vapour_mixing_ratio_and_theta_std.py's drho_dz,
# PySDM_examples/Shipway_and_Hill_2012/settings.py's drhod_dz). Unlike dP_dz
# above (which integrates P directly and needs no extra terms), PySDM integrates ρ_dry
# treating θ_dry/r as locally constant, plus a correction for the vertical mixing-ratio gradient dr/dz 
# θl_func is θ_dry(z), same convention as dP_dz
"""
dρ_dry_dz_pysdm(ρ_dry_z, ode_settings, z) 
solver for initial profiles computes the vertical dry air density gradient 
dρ_dry/dz at a given height z, based on the current dry air density ρ_dry_z
"""
function dρ_dry_dz_pysdm(ρ_dry_z, ode_settings, z)
    FT = eltype(ρ_dry_z)
    constants, θl_func, qt_func = ode_settings

    θ_dry_z::FT = θl_func(z)
    q_z::FT = qt_func(z)
    r_z::FT = mixing_ratio(q_z)

    T::FT = calc_T(constants, θ_dry_z, ρ_dry_z, q_z)
    P_dry::FT = ρ_dry_z * constants.Rd * T
    P::FT = calc_P_from_P_dry(P_dry, q_z, constants)

    Rq::FT = (constants.Rv * r_z + constants.Rd) / (1 + r_z)
    cp::FT = (constants.Cp_vapor * r_z + constants.Cp_air) / (1 + r_z)
    ρ_total::FT = P / (Rq * T)

    h = FT(1.0)
    dr_dz::FT = (mixing_ratio(qt_func(z + h)) - mixing_ratio(qt_func(z - h))) / (2h)

    drho_total_dz::FT = constants.gconst / T * ρ_total * (Rq / cp - 1) / Rq

    return drho_total_dz / (1 + r_z) - ρ_dry_z * dr_dz / (1 + r_z)^2
end

"""
init_droplets_scm(dist, settings, spatial, qv_profile; deterministic_multiplicity=false) 
initializes the single column model droplets based on the given distribution, coagulation settings, spatial settings, and water vapor profile.
"""
function init_droplets_scm(dist, settings::coag_settings{FT},
    spatial::spatial_settings{FT}, qv_profile::AbstractVector;
    deterministic_multiplicity::Bool=false) where FT<:AbstractFloat

    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0

    percentile_limit = [0.001, 0.999]
    ξstart::Vector{Int} = fill(1, Ns)
    z_loc = fill(FT(0.0), Ns)
    rad = fill(FT(0.0), Ns)
    cell_id = fill(1, Ns)

    nz = spatial.Nz

    # Allocate superdroplets proportionally to qv, only a thin layer above inversion (qt < 0.008)
    # Top half of BL gets 3x the weight of the bottom half
    weights = ifelse.(qv_profile .< FT(0.008), FT(0), qv_profile)

    bl_top = something(findlast(!iszero, weights), 0)
    mid = div(bl_top, 2)
    weights[mid+1:bl_top] .*= 3
    # seed a thin layer above the inversion (entrainment zone) so some aerosol exists
    # there too, at the same weight as the sub-cloud-deck layer (pre-boost)
    entrainment_zone_depth = FT(300.0) # m
    n_cells_above_inv = ceil(Int, entrainment_zone_depth / spatial.z_grid_height)
    above_inv = (bl_top+1):min(bl_top+n_cells_above_inv, nz)
    weights[above_inv] .= weights[1]
    weights ./= sum(weights)
    if spatial.weighted_droplet_allocation
        Ns_per_cell = floor.(Int, weights .* Ns)
        remainder = Ns - sum(Ns_per_cell)
        top_cells = sortperm(weights .* Ns .- Ns_per_cell, rev=true)[1:remainder]
        Ns_per_cell[top_cells] .+= 1
    else
        Ns_per_cell = fill(div(Ns, nz), nz)
    end

    offset = 0
    for k in 1:nz
        n_k = Ns_per_cell[k]
        drop_idx = offset+1 : offset+n_k

        cdf_values = if deterministic_multiplicity
            n_k == 0 ? FT[percentile_limit[1]] :
                collect(range(FT(percentile_limit[1]), FT(percentile_limit[2]), length=2*n_k + 1))
        else
            sort(rand(Uniform(percentile_limit[1], percentile_limit[2]), 2*n_k + 1))
        end
        rad[drop_idx] .= quantile.(dist, cdf_values[2:2:end-1])
        cdf_values = cdf_values[1:2:end]
        cdf_vals = cdf_values[2:end] - cdf_values[1:end-1]

        multiplicities = cdf_vals .* (n0 * ΔV)
        ξstart[drop_idx] .= floor.(Int, multiplicities .+ 0.5)
        z_loc[drop_idx] .= (k-1)*spatial.z_grid_height .+ spatial.z_grid_height .* rand(n_k)
        cell_id[drop_idx] .= k

        offset += n_k
    end

    dry_r3 = rad.^3
    volumes = 4* pi / 3 .* dry_r3 #initialize as just aerosol volume.

    dropidx = collect(1:Ns)

    sort!(dropidx, by = i -> cell_id[i])
    grid_range = map(1:spatial.Nz) do g
        s = findfirst(i -> cell_id[i] == g, dropidx)
        s === nothing ? (1:0) : s:findlast(i -> cell_id[i] == g, dropidx)
    end
    w_prime = zeros(FT, Ns) 

    droplets = droplet_attributes_1d{FT}(ξstart, volumes,dry_r3,z_loc, cell_id,w_prime,grid_range, dropidx)
    return droplets
end