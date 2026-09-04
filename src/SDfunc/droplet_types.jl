export droplet_attributes, droplet_attributes_1d, droplets_attributes_2d, simple_droplet_attributes
"""
    abstract type droplet_attributes{FT<:AbstractFloat}
"""
abstract type droplet_attributes{FT<:AbstractFloat} end

"""
    simple_droplet_attributes{FT<:AbstractFloat}
    struct for bare minimum droplet attributes, including multiplicity (ξ) and volume (X)
"""
struct simple_droplet_attributes{FT<:AbstractFloat} <:droplet_attributes{FT}
    ξ::Vector{Int}
    X::Vector{FT}
end

"""
    droplet_attributes_1d{FT<:AbstractFloat}
    struct for 1D droplet attributes, including multiplicity (ξ), volume (X), 
    dry volume (dry_r3), vertical location (z_loc), cell ID (cell_id), 
    vertical velocity perturbation (w_prime), grid range (grid_range), and index mapping (I)
"""
struct droplet_attributes_1d{FT<:AbstractFloat} <:droplet_attributes{FT}
    ξ::Vector{Int}
    X::Vector{FT}
    dry_r3::Vector{FT}
    # κ::Vector{FT}
    z_loc::Vector{FT}
    cell_id::Vector{Int}
    w_prime::Vector{FT}
    grid_range::Vector{UnitRange{Int}}
    I::Vector{Int}
end

"""
    droplet_attributes_2d{FT<:AbstractFloat}
    struct for 2D droplet attributes, including multiplicity (ξ), volume (X), 
    dry volume (dry_r3), vertical location (z_loc_in_cell), horizontal location (x_loc_in_cell), 
    and cell ID (cell_id), which is 1D index of the 2D grid cell.
    No corresponding 2D support yet
"""
struct droplet_attributes_2d{FT<:AbstractFloat} <:droplet_attributes{FT}
    ξ::Vector{Int}
    X::Vector{FT}
    dry_r3::Vector{FT}
    z_loc_in_cell::Vector{FT}
    x_loc_in_cell::Vector{FT}
    cell_id::Vector{Int}
end




droplet_attributes{FT}(ξ::Vector{Int}, X::Vector{FT}) where {FT<:AbstractFloat} = simple_droplet_attributes{FT}(ξ, X)
# droplet_attributes_1d{FT}(ξ::Vector{Int}, X::Vector{FT}, dry_mass::Vector{FT}, z_loc::Vector{FT}) where {FT<:AbstractFloat} = droplet_attributes_1d{FT}(ξ, X, dry_mass, z_loc)
droplet_attributes_2d{FT}(ξ::Vector{Int}, X::Vector{FT}, dry_mass::Vector{FT}, z_loc::Vector{FT}, x_loc::Vector{FT}) where {FT<:AbstractFloat} = droplet_attributes_2d{FT}(ξ, X, dry_mass, z_loc, x_loc)


Base.broadcastable(x::droplet_attributes_1d) = Ref(x)
Base.broadcastable(x::droplet_attributes_2d) = Ref(x)
Base.broadcastable(x::simple_droplet_attributes) = Ref(x)


