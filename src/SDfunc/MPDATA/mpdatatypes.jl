export BoundaryCondition, Periodic, NoFlux,
       limit, flux, mpdata_settings, mpdata_tmp, mpdata_fields, mpdata_mulitple_fields

abstract type BoundaryCondition end
struct Periodic <: BoundaryCondition end
struct NoFlux <: BoundaryCondition end

struct mpdata_settings
    n_corr::Int
    grid::Tuple{Int, Int}
    horizontal_boundary_condition::Union{Periodic, NoFlux}
    vertical_boundary_condition::Union{Periodic, NoFlux}

    function mpdata_settings(grid::Tuple{Int,Int};
                             n_corr=2,
                             horizontal_boundary_condition::BoundaryCondition=Periodic(),
                             vertical_boundary_condition::BoundaryCondition=Periodic())
        new(n_corr, grid, horizontal_boundary_condition, vertical_boundary_condition)
    end
end

struct mpdata_tmp
    ϕ::Matrix{Float64}
    GCx_step::Matrix{Float64}
    GCy_step::Matrix{Float64}
    GCx_tmp::Matrix{Float64}
    GCy_tmp::Matrix{Float64}

    function mpdata_tmp(ϕ, GCx_step, GCy_step)
        new(zeros(size(ϕ)), zeros(size(GCx_step)), zeros(size(GCy_step)), 
            zeros(size(GCx_step)), zeros(size(GCy_step)))
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
