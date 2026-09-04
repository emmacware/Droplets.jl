#---------------------------------------------------------
# Data Structures and Coalescence Simulation Timestep
#---------------------------------------------------------

#types
export Serial, Parallel,Adaptive,none,coagulation_run, coagulation_run_spatial
#coalescence sim
export coalescence_timestep!

# fix and export: collision_efficiency
struct Serial end
struct Parallel end
struct Adaptive end
struct none end
abstract type coagulation_run{FT<:AbstractFloat} end

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
struct coagulation_run_0D{FT<:AbstractFloat} <:coagulation_run{FT}
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

"""
    coagulation_run_spatial{FT<:AbstractFloat}

    struct for allocation related to coalescence in a spatially resolved simulation.
"""
struct coagulation_run_spatial{FT<:AbstractFloat} <:coagulation_run{FT}

    # N_in_cell::Vector{Int}
    Ngrids::Int
    I::Vector{Int}
    scale::Vector{FT} 
    # pαdt::Vector{Float64}
    # ϕ::Vector{Float64}
    lowest_zero::Ref{Bool}
    deficit::Vector{Float64}
    pair_starts::Vector{Int}  
    collision_rate::Vector{FT}
    collision_rate_pair::Vector{FT}  
    # grid_range::Vector{UnitRange{Int}}

    function coagulation_run_spatial{FT}(GridCount::Int, System_Ns,droplets) where FT<:AbstractFloat
        # N_in_cell = zeros(Int, GridCount)
        I = collect(1:System_Ns)
        scale = zeros(FT, GridCount)
        # pαdt = zeros(FT, div(System_Ns, 2))
        ϕ = zeros(FT, div(System_Ns, 2))
        lowest_zero = Ref(false)
        # deficit = zeros(FT, div(System_Ns, 2))
        deficit = zeros(FT, div(GridCount, 2))
        pair_starts = sizehint!(Int[], div(System_Ns, 2))
        collision_rate = zeros(FT, GridCount)
        collision_rate_pair = zeros(FT, System_Ns)


        new{FT}(GridCount,droplets.I,scale,lowest_zero, deficit, pair_starts, collision_rate, collision_rate_pair)
    end
end

Base.broadcastable(x::coagulation_run_spatial) = Ref(x)

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
- `run::backend`: Serial or Parallel
- `scheme::schemetype`: adaptive timestepping or none
- `droplets::droplet_attributes`: The superdroplets. 
if droplet_attributes_1d, coag_data must be coagulation_run_spatial

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



function rebuild_coalescence_pairing!(droplets::droplet_attributes_1d, coag_data::coagulation_run_spatial)
    empty!(coag_data.pair_starts)
    for g in eachindex(droplets.grid_range)
        isempty(droplets.grid_range[g]) && continue
        Ns_c = length(droplets.grid_range[g])
        coag_data.scale[g] = Ns_c * (Ns_c - 1) / 2 / div(Ns_c, 2)
        append!(coag_data.pair_starts, (droplets.grid_range[g])[1:2:end-1])
    end
    return nothing
end

function reduce_collision_rate!(droplets::droplet_attributes_1d, coag_data::coagulation_run_spatial)
    for i in coag_data.pair_starts
        cell = droplets.cell_id[coag_data.I[i]]
        coag_data.collision_rate[cell] += coag_data.collision_rate_pair[i]
    end
    return nothing
end

function coalescence_reshuffle_and_step!(::Serial, droplets::droplet_attributes_1d, coag_data::coagulation_run_spatial, kernel::Function, settings::coag_settings)
    for g in eachindex(droplets.grid_range)
        isempty(droplets.grid_range[g]) && continue
        shuffle!(@view coag_data.I[droplets.grid_range[g]])
    end
    for k in eachindex(coag_data.pair_starts)
        sdm_step!(coag_data.pair_starts[k], droplets, coag_data, kernel, settings)
    end
    reduce_collision_rate!(droplets, coag_data)
    return nothing
end

function coalescence_reshuffle_and_step!(::Parallel, droplets::droplet_attributes_1d, coag_data::coagulation_run_spatial, kernel::Function, settings::coag_settings)
    for g in eachindex(droplets.grid_range)
        isempty(droplets.grid_range[g]) && continue
        shuffle!(@view coag_data.I[droplets.grid_range[g]])
    end

    Threads.@threads for k in eachindex(coag_data.pair_starts)
        sdm_step!(coag_data.pair_starts[k], droplets, coag_data, kernel, settings)
    end
    reduce_collision_rate!(droplets, coag_data)
    return nothing
end

function coalescence_timestep!(::DynON,run::Union{Serial, Parallel},scheme::none,droplets::droplet_attributes_1d{FT},
    coag_data::coagulation_run_spatial,settings::coag_settings{FT}) where FT<:AbstractFloat
    rebuild_coalescence_pairing!(droplets, coag_data)
    coalescence_reshuffle_and_step!(run, droplets, coag_data, settings.kernel, settings)
    return nothing
end

function coalescence_timestep!(::DynOFF,run::Union{Serial, Parallel},scheme::none,droplets::droplet_attributes_1d{FT},
    coag_data::coagulation_run_spatial,settings::coag_settings{FT}) where FT<:AbstractFloat
end

function coalescence_timestep!(::DynON,run::Union{Serial, Parallel},scheme::none,droplets::droplet_attributes_1d{FT},
    coag_data::coagulation_run_spatial,settings::coag_settings{FT}, n_substeps::Int) where FT<:AbstractFloat
    rebuild_coalescence_pairing!(droplets, coag_data)
    for _ in 1:n_substeps
        coalescence_reshuffle_and_step!(run, droplets, coag_data, settings.kernel, settings)
    end
    return nothing
end

function coalescence_timestep!(::DynOFF,run::Union{Serial, Parallel},scheme::none,droplets::droplet_attributes_1d{FT},
    coag_data::coagulation_run_spatial,settings::coag_settings{FT}, n_substeps::Int) where FT<:AbstractFloat
end