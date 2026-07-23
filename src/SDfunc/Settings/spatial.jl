export spatial_settings, spatial_settings_1d, spatial_settings_2d, diagnostic_settings, scm_settings
export Dynamic, DynON, DynOFF
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
    area_per_grid::FT = 100.0 #m^2, only used for 1d case to calculate the volume of each grid cell
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

abstract type Dynamic end
struct DynON <: Dynamic end
struct DynOFF <: Dynamic end



struct scm_settings{FT<:AbstractFloat, Thr, Sch,
        Turb<:Dynamic, Cond<:Dynamic, REM_T<:Dynamic, Mot<:Dynamic,
        SpinupSat<:Dynamic, Rad<:Dynamic, Coag<:Dynamic, Sett<:Dynamic,
        Adv<:Dynamic, Rec<:Dynamic, TopEsc<:Dynamic, ThermoFB<:Dynamic, TurbDiff<:Dynamic}
    init_random_seed::Int
    coag_threading::Thr
    scheme::Sch
    Δt::FT
    n_cond::Int
    n_coag::Int
    spinup_time::FT
    turbulence::Turb
    condensation::Cond
    REM::REM_T
    motion::Mot
    spinupsaturation::SpinupSat
    radiation::Rad
    coalescence::Coag
    settling::Sett
    advection::Adv
    recycling::Rec
    top_escape::TopEsc
    thermo_feedback::ThermoFB
    turbulent_droplet_diffusion_on::TurbDiff
end

function scm_settings{FT}(;
        init_random_seed::Int             = Int(30),
        coag_threading                    = Serial(),
        scheme                            = none(),
        Δt::FT                            = FT(1.0),
        n_cond::Int                       = Int(10),
        n_coag::Int                       = Int(10),
        spinup_time::FT                   = FT(3600.0),
        turbulence::Dynamic               = DynON(),
        condensation::Dynamic             = DynON(),
        REM::Dynamic                      = DynON(),
        motion::Dynamic                   = DynON(),
        spinupsaturation::Dynamic         = DynON(),
        radiation::Dynamic                = DynON(),
        coalescence::Dynamic              = DynON(),
        settling::Dynamic                 = DynON(),
        advection::Dynamic                = DynON(),
        recycling::Dynamic                = DynON(),
        top_escape::Dynamic               = DynOFF(),
        thermo_feedback::Dynamic          = DynON(),
        turbulent_droplet_diffusion_on::Dynamic = DynON(),
    ) where {FT<:AbstractFloat}
    scm_settings{FT, typeof(coag_threading), typeof(scheme),
        typeof(turbulence), typeof(condensation), typeof(REM), typeof(motion),
        typeof(spinupsaturation), typeof(radiation), typeof(coalescence),
        typeof(settling), typeof(advection), typeof(recycling),
        typeof(top_escape), typeof(thermo_feedback), typeof(turbulent_droplet_diffusion_on)}(
        init_random_seed, coag_threading, scheme, Δt, n_cond, n_coag, spinup_time,
        turbulence, condensation, REM, motion, spinupsaturation, radiation,
        coalescence, settling, advection, recycling, top_escape, thermo_feedback,
        turbulent_droplet_diffusion_on)
end
