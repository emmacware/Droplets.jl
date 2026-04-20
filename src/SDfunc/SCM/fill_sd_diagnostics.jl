##Functions for turning superdroplet attributes into diagnostic variables in SCM
export sd_fill_diagnostics, scm_fill_diagnostic_output

function sd_fill_LWC(k::Int,sd::droplet_attributes{FT}, scm_grid,spatial::spatial_settings_1d,diagnostic_settings::diagnostic_settings) where {FT<:AbstractFloat}
    aerosol_grid_idxs = findall(sd.cell_id .== k .&& sd.X .< diagnostic_settings.aerosol_cloud_cuttoff)
    cloud_grid_idxs = findall(sd.cell_id .== k .&& sd.X .>= diagnostic_settings.aerosol_cloud_cuttoff .&& sd.X .< diagnostic_settings.cloud_rain_cuttoff)
    rain_grid_idxs = findall(sd.cell_id .== k .&& sd.X .>= diagnostic_settings.cloud_rain_cuttoff)

    if !isempty(aerosol_grid_idxs)
        scm_grid.diagnostics.aerosol_LWC[k] = sum(sd.X[aerosol_grid_idxs] .* sd.ξ[aerosol_grid_idxs]) .* constants.ρl
    end
    if !isempty(cloud_grid_idxs)
        scm_grid.diagnostics.cloud_LWC[k] = sum(sd.X[cloud_grid_idxs].* sd.ξ[cloud_grid_idxs]) .* constants.ρl
    end
    if !isempty(rain_grid_idxs)
        scm_grid.diagnostics.rain_LWC[k] = sum(sd.X[rain_grid_idxs] .* sd.ξ[rain_grid_idxs]) .* constants.ρl
    end 

    gridcell_volume = spatial.area_per_grid * spatial.z_grid_height
    scm_grid.diagnostics.aerosol_LWC[k] /= gridcell_volume
    scm_grid.diagnostics.cloud_LWC[k] /= gridcell_volume
    scm_grid.diagnostics.rain_LWC[k] /= gridcell_volume
end

function sd_fill_effective_radius(k::Int,sd::droplet_attributes{FT}, scm_grid,spatial::spatial_settings_1d,diagnostic_settings::diagnostic_settings) where {FT<:AbstractFloat}
    aerosol_grid_idxs = findall(sd.cell_id .== k .&& sd.X .< diagnostic_settings.aerosol_cloud_cuttoff)
    cloud_grid_idxs = findall(sd.cell_id .== k .&& sd.X .>= diagnostic_settings.aerosol_cloud_cuttoff .&& sd.X .< diagnostic_settings.cloud_rain_cuttoff)
    rain_grid_idxs = findall(sd.cell_id .== k .&& sd.X .>= diagnostic_settings.cloud_rain_cuttoff)

    if !isempty(aerosol_grid_idxs)
        aerosol_radii = volume_to_radius.(sd.X[aerosol_grid_idxs])
        scm_grid.diagnostics.aerosol_effective_radius[k] = sum(aerosol_radii.^3 .* sd.ξ[aerosol_grid_idxs]) / sum(aerosol_radii.^2 .* sd.ξ[aerosol_grid_idxs])
    end
    if !isempty(cloud_grid_idxs)
        cloud_radii = volume_to_radius.(sd.X[cloud_grid_idxs])
        scm_grid.diagnostics.cloud_effective_radius[k] = sum(cloud_radii.^3 .* sd.ξ[cloud_grid_idxs]) / sum(cloud_radii.^2 .* sd.ξ[cloud_grid_idxs])
    end
    if !isempty(rain_grid_idxs)
        rain_radii = volume_to_radius.(sd.X[rain_grid_idxs])
        scm_grid.diagnostics.rain_effective_radius[k] = sum(rain_radii.^3 .* sd.ξ[rain_grid_idxs]) / sum(rain_radii.^2 .* sd.ξ[rain_grid_idxs])
    end 

end

function sd_fill_diagnostics(sd::droplet_attributes{FT}, scm_grid,spatial::spatial_settings_1d,diagnostic_settings::diagnostic_settings) where {FT<:AbstractFloat}
    #set to zero before filling
        scm_grid.diagnostics.aerosol_LWC .= 0.0
        scm_grid.diagnostics.cloud_LWC .= 0.0
        scm_grid.diagnostics.rain_LWC .= 0.0
    
        scm_grid.diagnostics.aerosol_effective_radius .= 0.0
        scm_grid.diagnostics.cloud_effective_radius .= 0.0
        scm_grid.diagnostics.rain_effective_radius .= 0.0
    for k in 1:spatial.Nz
        sd_fill_LWC(k,sd, scm_grid,spatial, diagnostic_settings)
        sd_fill_effective_radius(k,sd, scm_grid,spatial, diagnostic_settings)
    end
end

function scm_fill_diagnostic_output(grid::scm_eulerian_arrays{FT},coagdata,conddata,raddata,spatial::spatial_settings_1d, t::Int) where FT
    grid.output.P[:,t] .= grid.states.P
    grid.output.qv[:,t] .= grid.states.qv
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
    # grid.output.surface_precipitation[t] done elsewhere

    return nothing
end

# lwc and eff_r for aerosol, cloud, rain