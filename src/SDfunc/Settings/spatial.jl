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
end

# Diagnostics Struct
Base.@kwdef struct diagnostic_settings{FT<:AbstractFloat}
    aerosol_cloud_cuttoff::FT = 1e-6
    cloud_rain_cuttoff::FT = 40e-6
end

Base.@kwdef struct scm_settings{FT<:AbstractFloat}

    init_random_seed::Int = Int(30) 
    coag_threading =  Serial()#Parallel(),use Julia NThreads for coalescence
    scheme = none() #Adaptive,Small_Alpha
    Δt::FT = FT(1.0)

end
