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
struct coag_settings{FT<:AbstractFloat, Kern}
    Δt::FT
    ΔV::FT
    Ns::Int
    scale::FT
    R_min::FT
    R_max::FT
    golovin_kernel_coeff::FT
    hydrodynamic_collision_eff_func::Bool
    long_lin_coeff::FT
    long_sq_coeff::FT
    long_r_thresh::FT
    kernel::Kern
    n0::FT
    R0::FT
end

function coag_settings{FT}(;
        Δt::FT                          = FT(1.0),
        ΔV::FT                          = FT(1e6),
        Ns::Int                         = 2^15,
        scale::FT                       = FT(Ns * (Ns - 1) / 2 / (Ns / 2)),
        R_min::FT                       = FT(1e-9),
        R_max::FT                       = FT(1e-3),
        golovin_kernel_coeff::FT        = FT(1.5e3),
        hydrodynamic_collision_eff_func::Bool = false,
        long_lin_coeff::FT              = FT(5.78e3),
        long_sq_coeff::FT               = FT(9.44e15),
        long_r_thresh::FT               = FT(5e-5),
        kernel                          = golovin,
        n0::FT                          = FT(2^23),
        R0::FT                          = FT(30.531e-6),
    ) where {FT<:AbstractFloat}
    coag_settings{FT, typeof(kernel)}(
        Δt, ΔV, Ns, scale, R_min, R_max, golovin_kernel_coeff,
        hydrodynamic_collision_eff_func, long_lin_coeff, long_sq_coeff,
        long_r_thresh, kernel, n0, R0)
end

Base.@kwdef struct condensation_settings{FT<:AbstractFloat}
    Δt::FT = FT(1.0) 
    activated::Bool = true # kohler or activated kohler, doesn't do anything yet
    #ammonium sulfate, Dziekan et al 2019 (from Petters and Kreidenweis (2007))
    kappa::FT = FT(0.61) # default Dycoms
    ρ_solute::FT = FT(1.78e3) # default Dycoms

end

Base.broadcastable(x::coag_settings) = Ref(x)
Base.broadcastable(x::condensation_settings) = Ref(x)