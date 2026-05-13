export spatial_settings, spatial_settings_1d, spatial_settings_2d, diagnostic_settings, scm_settings

abstract type spatial_settings{FT<:AbstractFloat} end

# Spatial Settings Struct
Base.@kwdef struct spatial_settings_2d{FT<:AbstractFloat} <:spatial_settings{FT}
    Nx::Int = 30 
    Nz::Int = 30 
    num_grids ::Int = Nx * Nz
    x_domain::FT = FT(1500.0) 
    z_domain::FT = FT(1500.0)
    z_grid_height::FT = z_domain / Nz 
    x_grid_width::FT = x_domain / Nx 
    periodic_boundaries_x::Bool = true 
    settling::Bool = true 
end

# Spatial Settings Struct
Base.@kwdef struct spatial_settings_1d{FT<:AbstractFloat} <:spatial_settings{FT}
    Nz::Int = 30 
    Z_max::FT = FT(1500.0)
    z_grid_height::FT = Z_max / Nz
    area_per_grid::FT = 1.0 #m^2, only used for 1d case to calculate the volume of each grid cell
    periodic_boundaries_x::Bool = true 
    settling::Bool = true 
    dt::FT = FT(1.0)
    t_max::Int = 3600
    dt_output::FT = FT(10.0)
end

# Diagnostics Struct
Base.@kwdef struct diagnostic_settings{FT<:AbstractFloat}
    aerosol_cloud_cuttoff::FT = radius_to_volume(1e-6)
    cloud_rain_cuttoff::FT = radius_to_volume(40e-6)
end

Base.@kwdef struct scm_settings{FT<:AbstractFloat}
    init_random_seed::Int = Int(30)
    coag_threading =  Serial()
    scheme = none()
    Δt::FT = FT(1.0)
    # surface_latent_heat_flux::FT = FT(93.0) # W/m^2
    # surface_sensible_heat_flux::FT = FT(16.0) # W/m^2
    n_cond::FT = Int(10)
    n_coag::FT = Int(10)
    spinup_time::FT = FT(3600.0)
    # Process on/off switches
    turbulence_on::Bool = true
    condensation_on::Bool = true
    radiation_on::Bool = true
    coalescence_on::Bool = true
    settling::Bool = true
    # Turbulence diffusion options
    # rho_weighted_diffusion::Bool = false          # use ρ-weighted diffuse_ρ_fields!; false = diffuse_fields!
    turbulent_droplet_diffusion_on::Bool = true   # run OU-process droplet w' kick
end
