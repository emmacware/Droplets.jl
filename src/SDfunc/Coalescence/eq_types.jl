#---------------------------------------------------------
# Data Structures and Coalescence Simulation Timestep
#---------------------------------------------------------

#types
export Serial, Parallel,Adaptive,none,coag_settings,coagulation_run, coagulation_run_spatial
#coalescence sim
export coalescence_timestep!

# fix and export: collision_efficiency
struct Serial end
struct Parallel end
struct Adaptive end
struct none end

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


"""
    struct coagulation_run{FT<:AbstractFloat}

Struct initialising temp memory used for coagulation.

Input:
- Ns::Int : Number of superdroplets.

Fields:
- I::Vector{Int} : Vector of indexes to be shuffled for permutations.
- pαdt::Vector{FT} : Vector of coalescence probabilities for each pair.
- ϕ::Vector{FT} : Random numbers to be used in Monte Carlo coalescence.
- lowest_zero::Ref{Bool} : Reference to a boolean indicating if the lowest multiplicity is zero.
- deficit::Ref{FT} : Reference to a float representing the deficit in mass or volume.

"""
struct coagulation_run{FT<:AbstractFloat}
    Ns::Int
    I::Vector{Int}
    scale::FT
    pαdt::Vector{FT}
    ϕ::Vector{FT}
    lowest_zero::Ref{Bool}
    deficit::Ref{FT}

    function coagulation_run{FT}(Ns::Int) where FT<:AbstractFloat
        I = collect(1:Ns)
        scale = div(Ns * (Ns - 1) , 2) / div(Ns , 2)
        pαdt = zeros(FT, div(Ns, 2))
        ϕ = zeros(FT, div(Ns, 2))
        lowest_zero = Ref(false)
        deficit = Ref(zero(FT))
        new{FT}(Ns,I, scale, pαdt, ϕ, lowest_zero, deficit)
    end
end

struct coagulation_run_spatial{FT<:AbstractFloat}

    N_in_cell::Vector{Int}
    I::Vector{Int}
    scale::Vector{FT}
    pαdt::Vector{Float64}
    ϕ::Vector{Float64}
    lowest_zero::Ref{Bool}
    deficit::Vector{Float64}

    function coagulation_run_spatial{FT}(GridCount::Int, System_Ns) where FT<:AbstractFloat
        N_in_cell = zeros(Int, div(System_Ns, GridCount))
        I = collect(1:System_Ns)
        pαdt = zeros(FT, div(System_Ns, 2))
        ϕ = zeros(FT, div(System_Ns, 2))
        lowest_zero = Ref(false)
        deficit = zeros(FT, length(N_in_cell))
        new{FT}(N_in_cell, I, pαdt, ϕ, lowest_zero, deficit)
    end
end



#----------------------------------------------------------
# COALESCENCE
#----------------------------------------------------------

"""
    coalescence_timestep!(run::backend, scheme::scheme_type, droplets::droplet_attributes)

Perform a coalescence timestep for the given droplets using the Superdroplet Method (SDM) 
Shima et al. (2009)
when the lowest multiplicity of superdroplets is less than 1, the 
largest superdroplet is split into two equal parts, as proposed by
Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017

# Arguments
- `run::backend`: Threading over Linear Sampling option
- `scheme::schemetype`: adaptive or none
- `droplets::droplet_attributes`: The superdroplets.

"""
function coalescence_timestep!(run::Union{Serial, Parallel},scheme::none, droplets::droplet_attributes,
    coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns::Int = settings.Ns
    
    shuffle!(coag_data.I)
    L = [(coag_data.I[l-1], coag_data.I[l]) for l in 2:2:Ns]

    compute_pαdt!(L, droplets,coag_data,settings.kernel,settings)

    rand!(coag_data.ϕ)

    test_pairs!(run,L,droplets,coag_data)

end 


function coalescence_timestep!(run::Union{Serial, Parallel},scheme::Adaptive,droplets::droplet_attributes{FT},
    coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    t_left = Ref(settings.Δt)

    while t_left[] > 0
        shuffle!(coag_data.I)
        L = [(coag_data.I[l-1], coag_data.I[l]) for l in 2:2:Ns]
        rand!(coag_data.ϕ)

        adaptive_pαdt!(L,droplets,coag_data,t_left,settings.kernel,settings)

        test_pairs!(run,L,droplets,coag_data)

    end
    return nothing
end 