export scm_states, scm_wind, scm_diagnostics, scm_eulerian_arrays, create_scm_grids

struct scm_states{FT<:AbstractFloat} #<:droplet_attributes{FT}
    P::Vector{FT}
    T::Vector{FT}
    θ::Vector{FT}
    qv::Vector{FT}
    ρ::Vector{FT}
    e::Vector{FT}   # turbulent kinetic energy (m²/s²)
end

struct scm_wind{FT<:AbstractFloat}
    u::Vector{FT}
    v::Vector{FT}
    w::Vector{FT}
end

struct scm_diagnostics{FT<:AbstractFloat}
    # rh::Vector{FT}
    aerosol_effective_radius::Vector{FT}
    cloud_effective_radius::Vector{FT}
    rain_effective_radius::Vector{FT}
    aerosol_LWC::Vector{FT}
    cloud_LWC::Vector{FT}
    rain_LWC::Vector{FT}
end

struct scm_eulerian_arrays{FT<:AbstractFloat}
    nz::Int
    dz::FT
    faces_z::Vector{FT}
    centers_z::Vector{FT}
    states::scm_states{FT}
    wind::scm_wind{FT}
    diagnostics::scm_diagnostics{FT}
end

function create_scm_grids(num_levels::Int, dz::FT) where {FT<:AbstractFloat}
    faces_z = collect(0:dz:(num_levels)*dz)
    centers_z = collect(dz/2:dz:((num_levels-1)*dz + dz/2))
    states = scm_states{FT}(zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels), fill(FT(1e-4), num_levels))
    wind = scm_wind{FT}(zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels+1))
    diagnostics = scm_diagnostics{FT}(zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels),
        zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels))
    return scm_eulerian_arrays{FT}(num_levels,dz,faces_z, centers_z, states, wind, diagnostics)
end