##Functions for turning superdroplet attributes into diagnostic variables in SCM
export sd_fill_diagnostics

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

# lwc and eff_r for aerosol, cloud, rain