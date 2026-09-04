export BoundaryCondition, Periodic, NoFlux, Reservoir, Extrapolated,
       limit, flux, mpdata_settings, mpdata_tmp, mpdata_fields, mpdata_mulitple_fields,
         mpdata_settings_1d, mpdata_tmp_1d, mpdata_fields_1d

"""
    BoundaryCondition
    Abstract type for boundary conditions in MPDATA advection scheme.
"""
abstract type BoundaryCondition end
struct Periodic <: BoundaryCondition end
struct NoFlux <: BoundaryCondition end
struct Reservoir <: BoundaryCondition end
struct Extrapolated <: BoundaryCondition end

"""
mpdata_settings
Struct for 2D MPDATA advection scheme Settings
    - n_corr::Int: Number of corrective steps in MPDATA
    - grid::Tuple{Int, Int}: Grid dimensions (Nx, Nz)
    - horizontal_boundary_condition::Union{Periodic, NoFlux, Extrapolated}: Horizontal boundary condition type
    - vertical_boundary_condition::Union{Periodic, NoFlux, Extrapolated}: Vertical boundary condition type
    - nonoscillatory::Bool: Flag for non-oscillatory option
    - infinite_gauge::Bool: Flag for infinite gauge option

"""
struct mpdata_settings
    n_corr::Int
    grid::Tuple{Int, Int}
    horizontal_boundary_condition::Union{Periodic, NoFlux, Extrapolated}
    vertical_boundary_condition::Union{Periodic, NoFlux, Extrapolated}
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

"""
mpdata_settings_1d{TVar<:ThermoVariable}
Struct for 1D MPDATA advection scheme Settings
    - n_corr::Int: Number of corrective steps in MPDATA
    - grid::Int: Grid dimension (Nz)
    - vertical_boundary_condition::Union{Periodic, NoFlux, Extrapolated}: Vertical boundary condition type
    - nonoscillatory::Bool: Flag for non-oscillatory option
    - infinite_gauge::Bool: Flag for infinite gauge option
    - thermo_variable::TVar: Type of thermodynamic variable (θ, qv)
    - topcell_boundary_condition::Union{Periodic, NoFlux, Reservoir}: Top cell boundary condition type (commented out)

"""
struct mpdata_settings_1d{TVar<:ThermoVariable}
    n_corr::Int
    grid::Int
    vertical_boundary_condition::Union{Periodic, NoFlux, Extrapolated}
    # topcell_boundary_condition::Union{Periodic, NoFlux, Reservoir}
    nonoscillatory::Bool
    infinite_gauge::Bool
    thermo_variable::TVar  # advect (θ,qv) directly, or (θ_l,qt) with θ/qv reconstructed after

    function mpdata_settings_1d(grid::Int;
                             n_corr=2,
                             vertical_boundary_condition::BoundaryCondition=Periodic(),
                             nonoscillatory::Bool=false,
                             infinite_gauge::Bool=false,
                             thermo_variable::ThermoVariable=ThetaQvVar())
        new{typeof(thermo_variable)}(n_corr, grid, vertical_boundary_condition,
            nonoscillatory, infinite_gauge, thermo_variable)
    end
end

"""
mpdata_tmp
mpdata_1d
Struct for MPData Allocations 
    - ϕ::Matrix{Float64}: Field variable
    - GCx_step::Matrix{Float64}: Step size in x-direction
    - GCy_step::Matrix{Float64}: Step size in y-direction
    - GCx_tmp::Matrix{Float64}: Temporary storage for x-direction
    - GCy_tmp::Matrix{Float64}: Temporary storage for y-direction
    - minmax::NamedTuple: Named tuple containing local minimum and maximum values
"""
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
    minmax::@NamedTuple{localmin::Vector{Float64}, localmax::Vector{Float64}}

    function mpdata_tmp_1d(ϕ, GCz_step)
        minmax = (localmin=zeros(length(ϕ)), localmax=zeros(length(ϕ)))
        new(zeros(length(ϕ)), zeros(length(GCz_step)), zeros(length(GCz_step)), minmax)
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