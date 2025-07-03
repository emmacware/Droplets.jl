
#---------------------------------------------------------
# Superdroplet Method logic
#---------------------------------------------------------
# SDM logic
export adaptive_pαdt!, pair_Ps_adaptive!,compute_pαdt!,split_highest_multiplicity!,test_pairs!,pair_Ps!,sdm_update!



"""
    pair_Ps_adaptive!(α::Int, pair::Tuple{Int,Int}, droplets::droplet_attributes, coag_data::coagulation_run, t_max::Vector{FT})

This function calculates the coalescence probability for a given pair of droplets, using adaptive timestepping logic (Bartman et al. 2021). 
It takes the following arguments:

- `α::Int`: The index of the coalescence model.
- `pair::Tuple{Int,Int}`: The indices of the droplets in the pair.
- `droplets::droplet_attributes`: The attributes of the droplets.
- `coag_data::coagulation_run`: The coagulation run data.
- `t_max::Vector{FT}`: The maximum timestep allowed given multiple sampling.

"""
function pair_Ps_adaptive!(α::Int,pair::Tuple{Int,Int}, droplets::droplet_attributes, coag_data::coagulation_run,t_max::Vector{FT},
    t_left::Ref{FT},kernel::Function,settings::coag_settings{FT}) where FT<:AbstractFloat
    if droplets.ξ[pair[1]] < droplets.ξ[pair[2]]
        k = pair[1]
        j = pair[2]
    else
        k = pair[2]
        j = pair[1]
    end

    ξj, ξk = droplets.ξ[j], droplets.ξ[k]

    tmp = ξj * kernel(droplets, pair, settings) * coag_data.scale / settings.ΔV

    coag_data.pαdt[α] = tmp
    if tmp*t_left[] > coag_data.ϕ[α]
        t_max[α] = (div(ξj, ξk)) / tmp # div gets the floor
    else
        t_max[α] = t_left[]
    end

end


"""
    compute_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes, coag_data::coagulation_run, kernel::Function, scale::FT,scalecoagsettings::coag_settings{FT}) where FT<:AbstractFloat

Map the probability function over the list of droplet pairs, L, and update the coagulation data in place.

# Arguments
- `L::Vector{Tuple{Int,Int}}`: List of droplet indices to be considered for coalescence.
- `droplets::droplet_attributes`: Droplet attributes.
- `coag_data::coagulation_run`: Coagulation data.
- `kernel::Function`: Coalescence kernel function.
- `coagsettings::coag_settings{FT}`: Coagulation settings.

"""
@inline function compute_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes,coag_data::coagulation_run,kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    map(i -> pair_Ps!(i, L[i], droplets,coag_data,kernel, coagsettings), eachindex(L))
    coag_data.pαdt .*= coag_data.scale * coagsettings.Δt / coagsettings.ΔV
end

@inline function pair_Ps!(α::Int, (j,k)::Tuple{Int,Int}, droplets::droplet_attributes,coag_data::coagulation_run,kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    coag_data.pαdt[α] = max(droplets.ξ[j], droplets.ξ[k]) * kernel(droplets,(j,k), coagsettings)
end


"""
    adaptive_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes, coag_data::coagulation_run, t_left::Ref{FT}, kernel::Function, settings::coag_settings{FT}) where FT<:AbstractFloat

Map the probability function over the list of droplet pairs, L, using adaptive timestepping logic and update the coagulation data in place.

# Arguments
- `L::Vector{Tuple{Int,Int}}`: List of droplet indices.
- `droplets::droplet_attributes`: Droplet attributes.
- `coag_data::coagulation_run`: Coagulation data.
- `t_left::Ref{FT}`: Timestep after adaptive substep.
- `kernel::Function`: Coalescence kernel function.
- `settings::coag_settings{FT}`: Coagulation settings.

"""
function adaptive_pαdt!(L::Vector{Tuple{Int,Int}},droplets::droplet_attributes,coag_data::coagulation_run,t_left::Ref{FT},kernel::Function, settings::coag_settings{FT}) where FT<:AbstractFloat
    t_max = Vector{FT}(undef, length(L))
    map(α -> pair_Ps_adaptive!(α,L[α],droplets,coag_data,t_max,t_left,kernel,settings),eachindex(coag_data.pαdt))
    Δtm = minimum(t_max)
    Δt = min(Δtm,t_left[])
    coag_data.pαdt .*= Δt
    t_left[] -= Δt
    return nothing
end



#----------------------------------------------------------
# UPDATE
#----------------------------------------------------------
"""
    test_pairs!(scheme::backend, Ns::Int, L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes, coag_data::coagulation_run)

Perform the SDM coalescence update for the superdroplets. Update the droplet attributes in place in case of 
    coalescence event.

# Arguments
- `scheme::Union{Serial, Parallel}`: The scheme to use for the update.
- `Ns::Int`: The number of superdroplets.
- `L::Vector{Tuple{Int,Int}}`: The list of droplet pairs.
- `droplets::droplet_attributes`: The droplet attributes.
- `coag_data::coagulation_run`: The coagulation data.
"""

function test_pairs!(scheme::Serial,L::Vector{Tuple{Int,Int}},droplets::droplet_attributes{FT},coag_data::coagulation_run) where FT<:AbstractFloat
    
    coag_data.lowest_zero[] = false
    for α::Int in 1:eachindex(L)
            
        if coag_data.ϕ[α] >= coag_data.pαdt[α]
            continue
        end
        sdm_update!(L[α], α, droplets,coag_data)
    end
    if coag_data.lowest_zero[] == true
        split_highest_multiplicity!(droplets)
    end
end

function test_pairs!(scheme::Parallel,L::Vector{Tuple{Int,Int}},droplets::droplet_attributes{FT},coag_data::coagulation_run) where FT<:AbstractFloat
    
    coag_data.lowest_zero[] = false
    Threads.@threads for α in 1:eachindex(L)
        if coag_data.ϕ[α] >= coag_data.pαdt[α]
            continue
        end
        sdm_update!(L[α], α, droplets,coag_data)
    end
    if coag_data.lowest_zero[] == true
        split_highest_multiplicity!(droplets)
    end
end

"""
    sdm_update!(pair::Tuple{Int,Int}, α::Int, droplets::droplet_attributes, coag_data::coagulation_run)
Perform the SDM coalescence update for the superdroplets when a coalesence event is determined, updating
    the droplet attributes in place.

# Arguments
- `pair::Tuple{Int,Int}`: The pair of droplets to be updated.
- `α::Int`: The index of the coalescence model.
- `droplets::droplet_attributes`: The droplet attributes.
- `coag_data::coagulation_run`: The coagulation data.
"""

function sdm_update!(pair::Tuple{Int,Int},α::Int, droplets::droplet_attributes{FT},coag_data::coagulation_run) where FT<:AbstractFloat

    j,k = pair
    if droplets.ξ[j] < droplets.ξ[k]
        j,k = k,j
    end

    ξj, ξk = droplets.ξ[j], droplets.ξ[k]
    pα =  coag_data.pαdt[α]

    pα_floor::FT = @fastmath floor(pα)
    γ::FT  = coag_data.ϕ[α] < pα - pα_floor ? pα_floor +1 : pα_floor

    if γ >= (floor_ξj_div_ξk = floor(ξj / ξk))
        coag_data.deficit[] += (γ - floor_ξj_div_ξk) * ξk
        γ = floor_ξj_div_ξk
    end

    if ξj > γ * ξk
        droplets.ξ[j] -= γ * ξk
        droplets.X[k] = γ * droplets.X[j] + droplets.X[k]
    else
        droplets.ξ[j] = floor(ξk / 2)
        droplets.ξ[k] -= droplets.ξ[j]
        droplets.X[k] = droplets.X[j] = γ *droplets.X[j] + droplets.X[k]

        if droplets.ξ[j] == 0
            coag_data.lowest_zero[] = true
        end
    end
    return nothing
end

"""
    split_highest_multiplicity!(droplets::droplet_attributes{FT}) where FT<:AbstractFloat
Split the superdroplet with the highest multiplicity into two equal parts, as proposed by
    Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
    This function is called when the lowest multiplicity of superdroplets is less than 1.
# Arguments
- `droplets::droplet_attributes{FT}`: The droplet attributes.
"""


function split_highest_multiplicity!(droplets::droplet_attributes{FT}) where FT<:AbstractFloat
    if maximum(droplets.ξ) > 1
        while (minimum(droplets.ξ) <= 0 && maximum(droplets.ξ) > 1)
            argmin_i, argmax_i = argmin(droplets.ξ), argmax(droplets.ξ)
            droplets.ξ[argmin_i] = floor(droplets.ξ[argmax_i]/2)
            droplets.X[argmin_i] = droplets.X[argmax_i]

            droplets.ξ[argmax_i] -= floor(droplets.ξ[argmax_i]/2)
        end
    elseif (maximum(droplets.ξ) <= 1)

        println("Highest superdroplet cannot be split")
        #right now, break the model until this situation gets handled
        if (maximum(droplets.ξ) < 1)
            error("Highest and Lowest Superdroplet have ξ==0")
        end

        #Later:remove superdroplet.. how to handle between cells?

        # # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
        # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)], ", removing from system")
        # deleteat!(R,argmin(ξ))
        # # deleteat!(M,argmin(ξ))
        # deleteat!(X,argmin(ξ))
        # deleteat!(ξ,argmin(ξ))
        # # Ns=Ns-1
    end
    return nothing
end




