##Functions for turning superdroplet attributes into diagnostic variables in SCM
export sd_fill_diagnostics, scm_fill_diagnostic_output

function sd_fill_diagnostics(sd::droplet_attributes{FT}, scm_grid::scm_eulerian_arrays{FT}, spatial::spatial_settings_1d, diagnostic_settings::diagnostic_settings) where {FT<:AbstractFloat}
    diag = scm_grid.diagnostics
    aero_cut = diagnostic_settings.aerosol_cloud_cuttoff
    cloud_cut = diagnostic_settings.cloud_rain_cuttoff
    scale = constants.ρl / (spatial.area_per_grid * spatial.z_grid_height)

    diag.aerosol_LWC .= 0
    diag.cloud_LWC .= 0
    diag.rain_LWC .= 0
    diag.aerosol_effective_radius .= 0
    diag.cloud_effective_radius .= 0
    diag.rain_effective_radius .= 0

    @inbounds for k in 1:spatial.Nz
        isempty(sd.grid_range[k]) && continue
        aero_Xξ = zero(FT); aero_r3ξ = zero(FT); aero_r2ξ = zero(FT)
        cloud_Xξ = zero(FT); cloud_r3ξ = zero(FT); cloud_r2ξ = zero(FT)
        rain_Xξ  = zero(FT); rain_r3ξ  = zero(FT); rain_r2ξ  = zero(FT)

        for idx in sd.grid_range[k]   # iterate range directly — no sd.I[range] temp array
            i   = sd.I[idx]
            X   = sd.X[i]
            ξ   = sd.ξ[i]
            r   = volume_to_radius(X)
            r2ξ = r * r * ξ
            r3ξ = r * r2ξ
            if X < aero_cut
                aero_Xξ += X * ξ; aero_r3ξ += r3ξ; aero_r2ξ += r2ξ
            elseif X < cloud_cut
                cloud_Xξ += X * ξ; cloud_r3ξ += r3ξ; cloud_r2ξ += r2ξ
            else
                rain_Xξ  += X * ξ; rain_r3ξ  += r3ξ; rain_r2ξ  += r2ξ
            end
        end

        diag.aerosol_LWC[k] = aero_Xξ * scale
        diag.cloud_LWC[k]   = cloud_Xξ * scale
        diag.rain_LWC[k]    = rain_Xξ  * scale
        aero_r2ξ  > 0 && (diag.aerosol_effective_radius[k] = aero_r3ξ  / aero_r2ξ)
        cloud_r2ξ > 0 && (diag.cloud_effective_radius[k]   = cloud_r3ξ / cloud_r2ξ)
        rain_r2ξ  > 0 && (diag.rain_effective_radius[k]    = rain_r3ξ  / rain_r2ξ)
    end
end

function scm_fill_diagnostic_output(grid::scm_eulerian_arrays{FT},coagdata::coagulation_run_spatial,conddata::condensation_data,raddata::Rad,spatial::spatial_settings_1d,constants::Constants{FT}, t::Int) where {FT,Rad}
    grid.output.P[:,t] .= grid.states.P
    grid.output.qv[:,t] .= grid.states.qv
    compute_ql_at_cell!.(grid.states, 1:grid.nz,constants)
    grid.output.ql[:,t] .= grid.states.ql_tmp
    grid.output.ρ[:,t] .= grid.states.ρ
    grid.output.θ[:,t] .= grid.states.θ
    grid.output.u[:,t] .= grid.wind.u
    grid.output.v[:,t] .= grid.wind.v
    # output.w[:,t] .= grid.states.w
    grid.output.e[:,t] .= grid.states.e
    grid.output.eps[:,t] .= grid.states.eps
    grid.output.aerosol_effective_radius[:,t] .= grid.diagnostics.aerosol_effective_radius
    grid.output.cloud_effective_radius[:,t] .= grid.diagnostics.cloud_effective_radius
    grid.output.rain_effective_radius[:,t] .= grid.diagnostics.rain_effective_radius
    grid.output.aerosol_LWC[:,t] .= grid.diagnostics.aerosol_LWC
    grid.output.cloud_LWC[:,t] .= grid.diagnostics.cloud_LWC
    grid.output.rain_LWC[:,t] .= grid.diagnostics.rain_LWC
    grid.output.collision_rate[:,t] .= coagdata.collision_rate / spatial.dt_output
    coagdata.collision_rate .= 0.0
    grid.output.condensation_src[:,t] .= conddata.condensation_src / spatial.dt_output
    conddata.condensation_src .= 0.0
    grid.output.condensation_rad_net[:,t] .= conddata.condensation_rad_net / spatial.dt_output
    conddata.condensation_rad_net .= 0.0
    grid.output.condensation_rad_abs[:,t] .= conddata.condensation_rad_abs / spatial.dt_output
    conddata.condensation_rad_abs .= 0.0
    grid.output.cloud_heating_rate[:,t] .= raddata.cloud_heating_delta / spatial.dt_output
    raddata.cloud_heating_delta .= 0.0
    grid.output.LWP[t] = sum(grid.diagnostics.cloud_LWC .+ grid.diagnostics.rain_LWC) * spatial.z_grid_height
    fillnum!(grid,t)
    # grid.output.surface_precipitation[t] done elsewhere

    return nothing
end

# lwc and eff_r for aerosol, cloud, rain

function fillnum!(grid::scm_eulerian_arrays{FT}, t::Int) where FT
    nz = grid.nz
    droplets = grid.states.droplets
    # X_threshold = FT(4π/3) * FT(1e-6)^3
    dziekan_cloud_fraction_bounds = (FT(5e-7), FT(25e-6))
    X_threshold_dziekan = radius_to_volume.(dziekan_cloud_fraction_bounds)
    for k in 1:nz
        r = droplets.grid_range[k]
        if isempty(r)
            grid.output.number[k,t] = 0.0
            continue
        end
        idxs = droplets.I[r]
        # grid.diagnostics.cloud_LWC[k] == 0 && continue


        # grid.output.number[k,t] = sum(droplets.ξ[i] for i in idxs if droplets.X[i] > X_threshold; init=0)
        # grid.output.number[k,t] = sum(droplets.ξ[idxs])# for i in idxs if droplets.X[i] > X_threshold; init=0)
        #now needs to be between dziekan bounds
        grid.output.number[k,t] = sum(droplets.ξ[i] for i in idxs if droplets.X[i] > X_threshold_dziekan[1]; init=0)
        grid.output.number[k,t] < 20*1e6*50 && (grid.output.number[k,t] = 0.0) 

    end
    return nothing
end