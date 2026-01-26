export spatial_settings

# Spatial Settings Struct
Base.@kwdef struct spatial_settings{FT<:AbstractFloat}
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