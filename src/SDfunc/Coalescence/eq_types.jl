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

struct coagulation_run_spatial{FT<:AbstractFloat} <:coagulation_run{FT}

    # N_in_cell::Vector{Int}
    Ngrids::Int
    I::Vector{Int}
    # scale::Vector{FT}
    # pαdt::Vector{Float64}
    # ϕ::Vector{Float64}
    lowest_zero::Ref{Bool}
    deficit::Vector{Float64}
    first_in_pair::Vector{Bool}
    # grid_range::Vector{UnitRange{Int}}

    function coagulation_run_spatial{FT}(GridCount::Int, System_Ns,droplets) where FT<:AbstractFloat
        # N_in_cell = zeros(Int, GridCount)
        I = collect(1:System_Ns)
        # scale = zeros(FT, div(System_Ns, 2))
        # pαdt = zeros(FT, div(System_Ns, 2))
        ϕ = zeros(FT, div(System_Ns, 2))
        lowest_zero = Ref(false)
        deficit = zeros(FT, div(System_Ns, 2))
        first_in_pair = falses(System_Ns)
        
        # grid_range = Vector{UnitRange{Int}}(undef, GridCount)
        # sort!(coagdata.I, by = i -> droplets.cell_id[i])
        # coagdata.grid_range .= [findfirst(i -> droplets.cell_id[i] == g, coagdata.I) : findlast(i -> droplets.cell_id[i] == g, coagdata.I) for g in 1:GridCount]

        new{FT}(GridCount,I,lowest_zero, deficit, first_in_pair)#, grid_range)
        # new{FT}(N_in_cell, I,scale, pαdt, ϕ, lowest_zero, deficit)
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


# #totally rework
# function coalescence_timestep!(run::Union{Serial, Parallel},scheme::none,droplets::droplet_attributes_1d{FT},
#     coag_data::coagulation_run_spatial,settings::coag_settings{FT}) where FT<:AbstractFloat
    
#     N_grids = length(coag_data.N_in_cell)
#     Ns = settings.Ns
#     shuffle!(coag_data.I)
#     #sort by cell id
#     sort!(coag_data.I, by = i -> droplets.cell_id[i])
#     #if there are d
#     coag_data.scale.= 0
#     #initiale L:
#     L = Vector{Tuple{Int,Int}}(undef, div(Ns, 2))

#     for g in 1:N_grids
#         #find the number of droplets in the grid
#         grid_idx = findfirst(i -> droplets.cell_id[i] == g, coag_data.I)
#         if grid_idx == nothing
#             continue
#         end
#         grid_Ns = count(i -> droplets.cell_id[i] == g, coag_data.I)
#         pair_idx_start = div(grid_idx, 2) + 1
#         pair_idx_end = pair_idx_start + div(grid_Ns-1, 2) - 1
        
#         #find first empty pair in L:
#         first_empty_pair_idx = findfirst(i -> L[i] == nothing, 1:length(L))
#         #fill L with pairs of droplets in the grid
#         for i in 0:div(grid_Ns-1, 2)-1
#             L[first_empty_pair_idx + i] = (coag_data.I[grid_idx + 2*i - 1], coag_data.I[grid_idx + 2*i])
#         end

#         # fill scale for index grid_idx/2 to grid_idx/2 + grid_num//2 with the appropriate value
#         coag_data.scale[first_empty_pair_idx:pair_idx_end] .= div(grid_Ns * (grid_Ns - 1) , 2) / div(grid_Ns , 2)

#     end


#     compute_pαdt!(L, droplets,coag_data,settings.kernel,settings) # check if this still works with the scale being a vector

#     rand!(coag_data.ϕ)

#     test_pairs!(run,L,droplets,coag_data)

#     return nothing
# end

#totally rework
function coalescence_timestep!(run::Union{Serial, Parallel},scheme::none,droplets::droplet_attributes_1d{FT},
    coag_data::coagulation_run_spatial,settings::coag_settings{FT}) where FT<:AbstractFloat
    
    # N_grids = coag_data.Ngrids
    Ns = settings.Ns
    coag_data.first_in_pair .= false
    sort!(coag_data.I, by = i -> droplets.cell_id[i])
    # droplets.grid_range .= [findfirst(i -> droplets.cell_id[i] == g, coag_data.I) : findlast(i -> droplets.cell_id[i] == g, coag_data.I) for g in 1:coag_data.Ngrids]

    for g in eachindex(droplets.grid_range)
        start = findfirst(i -> droplets.cell_id[i] == g, coag_data.I) 
        if start == nothing
            continue
        end
        droplets.grid_range[g] = start:findlast(i -> droplets.cell_id[i] == g, coag_data.I)
        coag_data.I[droplets.grid_range[g]] .= shuffle(coag_data.I[droplets.grid_range[g]]) 
        coag_data.first_in_pair[(droplets.grid_range[g])[1:2:end-1]] .= true

    end

    map(i -> sdm_step!(i,droplets,coag_data,settings.kernel, settings), 1:Ns)

    # compute_pαdt!(droplets,coag_data,settings.kernel,settings) # check if this still works with the scale being a vector

    # rand!(coag_data.ϕ)

    # test_pairs!(run,L,droplets,coag_data)

    return nothing
end