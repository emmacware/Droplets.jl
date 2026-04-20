export coag_settings, condensation_settings

"""
    coag_settings{FT<:AbstractFloat}

A struct representing the settings for coalescence.

# Fields
- `FT`: AbstractFloat.
- `Δt`: The time step for the coalescence simulation.
- `ΔV`: The volume of the coalescence domain.
- `Ns`: The number of superdroplets.
- `R_min`: The minimum radius of sampled droplets, used in some initializations. Radius in meters
- `R_max`: The maximum radius of sampled droplets, used in some initializations. Radius in meters.
- `golovin_kernel_coeff`: The Golovin kernel coefficient, 1/seconds.
- `hydrodynamic_collision_eff_func`: A boolean indicating whether to use the
        hydrodynamic collision efficiency function. Currently not implemented.
- `kernel`: The kernel function used for coalescence. Implemented options are `golovin` or `hydrodynamic`,
    can take any user implemented function with input type (::droplet_attributes, pairindex::Tuple{Int,Int}, ::coag_settings{FT}).
- `n0`: The initial real world droplet concentration.
- `R0`: The initial seed radius of the droplets.
"""
Base.@kwdef struct coag_settings{FT<:AbstractFloat}
    Δt::FT = FT(1.0)
    ΔV::FT = FT(1e6)
    Ns::Int = 2^15# number of superdroplets
    scale::FT = Ns * (Ns - 1) / 2 / (Ns / 2)
    R_min::FT = FT(1e-9)
    R_max::FT = FT(1e-3)
    golovin_kernel_coeff::FT = FT(1.5e3)
    hydrodynamic_collision_eff_func::Bool = false
    kernel::Function = golovin # golovin, hydrodynamic
    n0::FT = FT(2^23) # initial droplet concentration
    R0::FT = FT(30.531e-6) # initial radius
end

Base.@kwdef struct condensation_settings{FT<:AbstractFloat}
    Δt::FT = FT(1.0) 
    activated::Bool = true # kohler or activated kohler, doesn't do anything yet
    kappa::FT = FT(0.0) # hygroscopicity parameter
    ρ_solute::FT = FT(1.78e3) # default Dycoms

end