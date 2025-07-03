using StaticArrays

export droplet_attributes, droplet_attributes_1d, droplets_attributes_2d, simple_droplet_attributes, static_droplet_attributes
"""
    abstract type droplet_attributes{FT<:AbstractFloat}
"""
abstract type droplet_attributes{FT<:AbstractFloat} end


struct simple_droplet_attributes{FT<:AbstractFloat} <:droplet_attributes{FT}
    ξ::Vector{Int}
    X::Vector{FT}
end

struct static_droplet_attributes{FT<:AbstractFloat,NSD} <:droplet_attributes{FT}
    ξ::SVector{NSD,FT}
    X::SVector{NSD,FT}
end

struct droplet_attributes_1d{FT<:AbstractFloat} <:droplet_attributes{FT}
    ξ::Vector{Int}
    X::Vector{FT}
    dry_r3::Vector{FT}
    z_loc_in_cell::Vector{FT}
    cell_id::Vector{Int}
end

struct droplet_attributes_2d{FT<:AbstractFloat} <:droplet_attributes{FT}
    ξ::Vector{Int}
    X::Vector{FT}
    dry_r3::Vector{FT}
    z_loc_in_cell::Vector{FT}
    x_loc_in_cell::Vector{FT}
    cell_id::Vector{Int}
end




"""
droplet_attributes{FT} where {FT<:AbstractFloat}
Create a new instance of simple_droplet_attributes with the given attribute vectors.
"""
droplet_attributes{FT}(ξ::Vector{Int}, X::Vector{FT}) where {FT<:AbstractFloat} = simple_droplet_attributes{FT}(ξ, X)
droplet_attributes{FT}(ξ::SVector{NSD,FT}, X::SVector{NSD,FT}) where {FT<:AbstractFloat, NSD} = static_droplet_attributes{FT, NSD}(ξ, X)
droplet_attributes_1d{FT}(ξ::Vector{Int}, X::Vector{FT}, dry_mass::Vector{FT}, z_loc::Vector{FT}) where {FT<:AbstractFloat} = droplet_attributes_1d{FT}(ξ, X, dry_mass, z_loc)
droplet_attributes_2d{FT}(ξ::Vector{Int}, X::Vector{FT}, dry_mass::Vector{FT}, z_loc::Vector{FT}, y_loc::Vector{FT}) where {FT<:AbstractFloat} = droplet_attributes_2d{FT}(ξ, X, dry_mass, z_loc, y_loc)




