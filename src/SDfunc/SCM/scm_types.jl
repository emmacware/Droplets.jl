export scm_states, scm_wind, scm_diagnostics, scm_eulerian_arrays, create_scm_grids
export scm_outputs, condensation_data #, radiation_data
export turbulence_data

struct scm_states{FT<:AbstractFloat} #<:droplet_attributes{FT}
    nz::Int
    P::Vector{FT}
    P_faces::Vector{FT}
    T_tmp::Vector{FT}
    ql_tmp::Vector{FT}
    θ::Vector{FT}
    qv::Vector{FT}
    ρ::Vector{FT}
    e::Vector{FT}   # turbulent kinetic energy (m²/s²)
    eps::Vector{FT} # TKE dissipation rate (m²/s³)
    droplets::droplet_attributes_1d{FT}
    spatial::spatial_settings_1d{FT}
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

struct scm_outputs{FT<:AbstractFloat}
    time::Vector{FT}
    P::Matrix{FT}
    qv::Matrix{FT}
    ql::Matrix{FT}
    ρ::Matrix{FT}
    θ::Matrix{FT}
    u::Matrix{FT}
    v::Matrix{FT}
    # w::Matrix{FT}
    e::Matrix{FT}
    eps::Matrix{FT}
    aerosol_effective_radius::Matrix{FT}
    cloud_effective_radius::Matrix{FT}
    rain_effective_radius::Matrix{FT}
    aerosol_LWC::Matrix{FT}
    cloud_LWC::Matrix{FT}
    rain_LWC::Matrix{FT}
    collision_rate::Matrix{FT}
    condensation_src::Matrix{FT}
    condensation_rad_net::Matrix{FT}
    condensation_rad_abs::Matrix{FT}
    cloud_heating_rate::Matrix{FT}

    surface_precipitation::Vector{FT}
    LWP::Vector{FT}
end

function scm_outputs(num_levels::Int, t_max::Int, dt_output::FT)::scm_outputs{FT} where FT<:AbstractFloat
    time = collect(0:dt_output:t_max)
    n_output_steps = length(time)
    P = zeros(FT, num_levels, n_output_steps)
    qv = zeros(FT, num_levels, n_output_steps)
    ql = zeros(FT, num_levels, n_output_steps)
    ρ = zeros(FT, num_levels, n_output_steps)
    θ = zeros(FT, num_levels, n_output_steps)
    u = zeros(FT, num_levels, n_output_steps)
    v = zeros(FT, num_levels, n_output_steps)
    # w = zeros(FT, num_levels, n_output_steps)
    e = zeros(FT, num_levels, n_output_steps)
    eps = zeros(FT, num_levels, n_output_steps)
    aerosol_effective_radius = zeros(FT, num_levels, n_output_steps)
    cloud_effective_radius = zeros(FT, num_levels, n_output_steps)
    rain_effective_radius = zeros(FT, num_levels, n_output_steps)
    aerosol_LWC = zeros(FT, num_levels, n_output_steps)
    cloud_LWC = zeros(FT, num_levels, n_output_steps)
    rain_LWC = zeros(FT, num_levels, n_output_steps)
    collision_rate = zeros(FT, num_levels, n_output_steps)
    condensation_src = zeros(FT, num_levels, n_output_steps)
    condensation_rad_net = zeros(FT, num_levels, n_output_steps)
    condensation_rad_abs = zeros(FT, num_levels, n_output_steps)
    cloud_heating_rate = zeros(FT, num_levels, n_output_steps)   

    surface_precipitation = zeros(FT, n_output_steps)
    LWP = zeros(FT, n_output_steps)

    return scm_outputs{FT}(time,P,qv,ql,ρ,θ,u,v,e,eps,aerosol_effective_radius,
        cloud_effective_radius,rain_effective_radius,aerosol_LWC,
        cloud_LWC,rain_LWC,
        collision_rate,
        condensation_src,
        condensation_rad_net,
        condensation_rad_abs,
        cloud_heating_rate,
        surface_precipitation,
        LWP
        )
end

struct scm_eulerian_arrays{FT<:AbstractFloat}
    nz::Int
    dz::FT
    faces_z::Vector{FT}
    centers_z::Vector{FT}
    states::scm_states{FT}
    wind::scm_wind{FT}
    diagnostics::scm_diagnostics{FT}
    spatial::spatial_settings_1d{FT}
    output::scm_outputs{FT}
end

function create_scm_grids(num_levels::Int, dz::FT,droplets::droplet_attributes_1d;spatial::Union{spatial_settings_1d,Nothing}=nothing,output::Union{scm_outputs,Nothing}=nothing) where {FT<:AbstractFloat}
    faces_z = collect(0:dz:(num_levels)*dz)
    centers_z = collect(dz/2:dz:((num_levels-1)*dz + dz/2))
    wind = scm_wind{FT}(zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels+1))
    diagnostics = scm_diagnostics{FT}(zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels),
        zeros(FT,num_levels), zeros(FT,num_levels), zeros(FT,num_levels))
    if spatial == nothing
        spatial = spatial_settings_1d{FT}(Nz=num_levels, z_grid_height=dz, Z_max=num_levels*dz)
    end
    states = scm_states{FT}(num_levels,zeros(FT,num_levels),zeros(FT,num_levels+1), zeros(FT,num_levels), 
        zeros(FT,num_levels),zeros(FT,num_levels),zeros(FT,num_levels),zeros(FT,num_levels), zeros(FT,num_levels),zeros(FT,num_levels), droplets,spatial)
    if output == nothing
        output = scm_outputs(num_levels, spatial.t_max, spatial.dt_output)
    end
    return scm_eulerian_arrays{FT}(num_levels,dz,faces_z, centers_z, states, wind, diagnostics,spatial,output)
end


struct condensation_data{FT<:AbstractFloat}
    condensation_src::Vector{FT}
    condensation_rad_net::Vector{FT}
    condensation_rad_abs::Vector{FT}
    vol_change_helper::Vector{FT}
end

function condensation_data(::Type{FT},num_levels::Int,nsd::Int)::condensation_data{FT} where FT<:AbstractFloat
    condensation_src = zeros(FT, num_levels)
    condensation_rad_net = zeros(FT, num_levels)
    condensation_rad_abs = zeros(FT, num_levels)
    vol_change_helper = zeros(FT, num_levels)
    return condensation_data{FT}(condensation_src, condensation_rad_net, condensation_rad_abs,vol_change_helper)
end

struct turbulence_data{FT<:AbstractFloat} 
    l::Vector{FT} # mixing length
    SM::Vector{FT} # shear production
    SH::Vector{FT} # buoyancy production
    GM::Vector{FT} # shear diffusion
    GH::Vector{FT} # buoyancy diffusion
    K_h::Vector{FT} # diffusivity for scalars
    K_m::Vector{FT} # diffusivity for momentum
    K_e::Vector{FT} # diffusivity for TKE
    du::Vector{FT} # 
    dv::Vector{FT} #

    function turbulence_data(::Type{FT}, num_levels::Int)::turbulence_data{FT} where FT<:AbstractFloat
        l = zeros(FT, num_levels)
        SM = zeros(FT, num_levels)
        SH = zeros(FT, num_levels)
        GM = zeros(FT, num_levels)
        GH = zeros(FT, num_levels)
        K_h = zeros(FT, num_levels)
        K_m = zeros(FT, num_levels)
        K_e = zeros(FT, num_levels)
        du = zeros(FT, num_levels)
        dv = zeros(FT, num_levels)
        return new{FT}(l, SM, SH, GM, GH, K_h, K_m, K_e, du, dv)
    end
end


Base.broadcastable(x::scm_states) = Ref(x)
Base.broadcastable(x::scm_wind) = Ref(x)
Base.broadcastable(x::scm_diagnostics) = Ref(x)
Base.broadcastable(x::scm_eulerian_arrays) = Ref(x)
Base.broadcastable(x::scm_outputs) = Ref(x)