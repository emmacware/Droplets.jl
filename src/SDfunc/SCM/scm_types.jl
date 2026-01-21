export scm_states, scm_wind, scm_diagnostics, scm_eulerian_arrays, create_scm_grids

struct scm_states{FT<:AbstractFloat} #<:droplet_attributes{FT}
    P::Vector{FT}
    T::Vector{FT}
    θ::Vector{FT}
    qv::Vector{FT}
end

struct scm_wind{FT<:AbstractFloat}
    u::Vector{FT}
    v::Vector{FT}
    w::Vector{FT}
end

struct scm_diagnostics{FT<:AbstractFloat}
    # rh::Vector{FT}
    effective_radius::Vector{FT}
    nc::Vector{FT}
    LWC::Vector{FT}
end

struct scm_eulerian_arrays{FT<:AbstractFloat}
    faces_z::Vector{FT}
    centers_z::Vector{FT}
    states::scm_states{FT}
    wind::scm_wind{FT}
    diagnostics::scm_diagnostics{FT}
end

function create_scm_grids(num_levels::Int, dz::FT,halo::Int) where {FT<:AbstractFloat}
    faces_z = collect(0:dz:(num_levels)*dz)
    centers_z = collect(dz/2:dz:(num_levels-1)*dz + dz/2)
    states = scm_states{FT}(zeros(FT,num_levels+2*halo), zeros(FT,num_levels+2*halo), zeros(FT,num_levels+2*halo), zeros(FT,num_levels+2*halo))
    wind = scm_wind{FT}(zeros(FT,num_levels+2*halo), zeros(FT,num_levels+2*halo), zeros(FT,num_levels+1+2*halo))
    diagnostics = scm_diagnostics{FT}(zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels))
    return scm_eulerian_arrays{FT}(faces_z, centers_z, states, wind, diagnostics)
end