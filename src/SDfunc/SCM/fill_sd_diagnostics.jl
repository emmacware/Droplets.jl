##Functions for turning superdroplet attributes into diagnostic variables in SCM
export sd_fill_diagnostics, scm_fill_diagnostic_output

function sd_fill_diagnostics(sd::droplet_attributes{FT}, scm_grid::scm_eulerian_arrays{FT}, spatial::spatial_settings_1d, diagnostic_settings::diagnostic_settings) where {FT<:AbstractFloat}
    diag = scm_grid.diagnostics
    aero_cut = diagnostic_settings.aerosol_cloud_cuttoff
    cloud_cut = diagnostic_settings.cloud_rain_cuttoff
    scale = 1 / (spatial.area_per_grid * spatial.z_grid_height)
    massscale = constants.ρl * scale

    diag.aerosol_LWC .= 0
    diag.cloud_LWC .= 0
    diag.rain_LWC .= 0
    diag.aerosol_effective_radius .= 0
    diag.cloud_effective_radius .= 0
    diag.rain_effective_radius .= 0

    @inbounds for k in 1:spatial.Nz
        isempty(sd.grid_range[k]) && continue
        aero_ξ = zero(FT); cloud_ξ = zero(FT); rain_ξ = zero(FT)
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
                aero_Xξ += X * ξ; aero_r3ξ += r3ξ; aero_r2ξ += r2ξ; aero_ξ += ξ
                
            elseif X < cloud_cut
                cloud_Xξ += X * ξ; cloud_r3ξ += r3ξ; cloud_r2ξ += r2ξ; cloud_ξ += ξ
            else
                rain_Xξ  += X * ξ; rain_r3ξ  += r3ξ; rain_r2ξ  += r2ξ; rain_ξ  += ξ
            end
        end

        diag.aerosol_LWC[k] = aero_Xξ * massscale
        diag.cloud_LWC[k]   = cloud_Xξ * massscale
        diag.rain_LWC[k]    = rain_Xξ  * massscale
        diag.aerosol_number[k] = aero_ξ * scale
        diag.cloud_number[k]   = cloud_ξ * scale
        diag.rain_number[k]    = rain_ξ  * scale
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
    grid.output.aerosol_number[:,t] .= grid.diagnostics.aerosol_number
    grid.output.cloud_number[:,t] .= grid.diagnostics.cloud_number
    grid.output.rain_number[:,t] .= grid.diagnostics.rain_number
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
    # fillnum!(grid,t)
    # grid.output.surface_precipitation[t] done elsewhere

    return nothing
end
