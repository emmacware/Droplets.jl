export BoundaryCondition, Periodic, NoFlux,
       limit, flux, mpdata_settings, mpdata_tmp, mpdata_fields, mpdata_mulitple_fields,
         mpdata_settings_1d, mpdata_tmp_1d, mpdata_fields_1d


abstract type BoundaryCondition end
struct Periodic <: BoundaryCondition end
struct NoFlux <: BoundaryCondition end

struct mpdata_settings
    n_corr::Int
    grid::Tuple{Int, Int}
    horizontal_boundary_condition::Union{Periodic, NoFlux}
    vertical_boundary_condition::Union{Periodic, NoFlux}
    nonoscillatory::Bool
    infinite_gauge::Bool

    function mpdata_settings(grid::Tuple{Int,Int};
                             n_corr=2,
                             horizontal_boundary_condition::BoundaryCondition=Periodic(),
                             vertical_boundary_condition::BoundaryCondition=Periodic(),
                             nonoscillatory::Bool=false,
                             infinite_gauge::Bool=false)
        new(n_corr, grid, horizontal_boundary_condition, vertical_boundary_condition,
            nonoscillatory, infinite_gauge)
    end
end

struct mpdata_settings_1d
    n_corr::Int
    grid::Int
    vertical_boundary_condition::Union{Periodic, NoFlux}
    nonoscillatory::Bool
    infinite_gauge::Bool

    function mpdata_settings_1d(grid::Int;
                             n_corr=2,
                             vertical_boundary_condition::BoundaryCondition=Periodic(),
                             nonoscillatory::Bool=false,
                             infinite_gauge::Bool=false)
        new(n_corr, grid, vertical_boundary_condition,
            nonoscillatory, infinite_gauge)
    end
end

struct mpdata_tmp
    ϕ::Matrix{Float64}
    GCx_step::Matrix{Float64}
    GCy_step::Matrix{Float64}
    GCx_tmp::Matrix{Float64}
    GCy_tmp::Matrix{Float64}
    minmax::@NamedTuple{localmin::Matrix{Float64}, localmax::Matrix{Float64}}

    function mpdata_tmp(ϕ, GCx_step, GCy_step)
        minmax = (localmin=zeros(size(ϕ)), localmax=zeros(size(ϕ)))
        new(zeros(size(ϕ)), zeros(size(GCx_step)), zeros(size(GCy_step)), 
            zeros(size(GCx_step)), zeros(size(GCy_step)),minmax)
    end
end

struct mpdata_tmp_1d
    ϕ::Vector{Float64}
    GCz_step::Vector{Float64}
    GCz_tmp::Vector{Float64}
    # minmax::@NamedTuple{localmin::Vector{Float64}, localmax::Vector{Float64}}

    function mpdata_tmp_1d(ϕ, GCz_step)
        # minmax = (localmin=zeros(size(ϕ)), localmax=zeros(size(ϕ)))
        new(zeros(length(ϕ)), zeros(length(GCz_step)), zeros(length(GCz_step)), 
            # minmax
            )
    end
end

struct mpdata_fields
    ϕ::Matrix{Float64}
    GCx::Matrix{Float64}
    GCy::Matrix{Float64}
end
struct mpdata_mulitple_fields
    ϕ::Tuple{Matrix{Float64}}
    GCx::Matrix{Float64}
    GCy::Matrix{Float64}
end

struct mpdata_fields_1d
    ϕ::Vector{Float64}
    GCz::Vector{Float64}
end