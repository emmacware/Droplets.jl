using Distributions
using Droplets
using Plots

function dP_dz(P_z, ode_settings, z)

    FT = eltype(P_z)
    constants,θl_func,qt_func = ode_settings

    θ_z::FT = θl_func(z)
    q_vap::FT = qt_func(z)

    g::FT = constants.gconst
    R_d::FT = constants.Rd
    R_v::FT = constants.Rv
    cp_d::FT = constants.Cp_air
    cp_v::FT = constants.Cp_vapor #double check

    T::FT = T_from_theta(θ_z,P_z,constants)
    T_virtual = T * (1 + 0.61 * q_vap)
    ρ::FT = P_z / (R_d * T_virtual)

    return -g * ρ
end


function initialize_scm_environment(nz, dz, P_surface, θl, qt, prescribed_u, prescribed_v, prescribed_w)
    grid = create_scm_grids(nz, dz)
    p = (constants, θl, qt)

    press = ODEProblem(dP_dz,P_surface,(0,maximum(grid.centers_z)+1),p)
    P_init = solve(press,Tsit5(), reltol = 1e-10, abstol = 1e-10,saveat = grid.centers_z)
    grid.states.P .= P_init.u

    grid.states.θ .= θl.(grid.centers_z)
    grid.states.qv .= qt.(grid.centers_z)
    grid.states.T .= T_from_theta(grid.states.θ,grid.states.P,constants)
    # grid.states.ρ .= grid.states.P ./ (constants.Rd .* T_virtual.(grid.states.T, grid.states.qv))
    grid.states.ρ .= ρ_ideal_gas(grid.states.P, grid.states.T, grid.states.qv, constants)

    grid.wind.u .= prescribed_u.(grid.centers_z)
    grid.wind.v .= prescribed_v.(grid.centers_z)
    grid.wind.w .= prescribed_w.(grid.faces_z)
    return grid
end


function init_droplets_dycoms_scm(dist,settings::coag_settings{FT},
    spatial::spatial_settings{FT})where FT<:AbstractFloat

    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0


    percentile_limit = [0.001, 0.999]  
    Rmin = quantile(dist, percentile_limit[1])
    Rmax = quantile(dist, percentile_limit[2])   

    Rarray = range(Rmin, Rmax, length=Ns+1)
    cdf_values = cdf.(dist, Rarray)
    rad = (Rarray[2:end] .+ Rarray[1:end-1])./ 2 
    cdf_values = cdf_values[2:end] - cdf_values[1:end-1]
    multiplicities = cdf_values .* (n0 * spatial.Z_max * spatial.area_per_grid)
                    
    ξstart::Vector{Int} = floor.(multiplicities .+ 0.5)
    z_loc = rand(Uniform(0, 1), Ns) * spatial.Z_max
    cell_id = Int.(floor.(z_loc ./ spatial.z_grid_height) .+ 1)

    dry_r3 = rad.^3
    volumes = 4* pi / 3 .* dry_r3 #initialize as just aerosol volume.

    dropidx = collect(1:Ns)

    sort!(dropidx, by = i -> cell_id[i])
    grid_range = [findfirst(i -> cell_id[i] == g, dropidx) : findlast(i -> cell_id[i] == g, dropidx) for g in 1:spatial.Nz]


    droplets = droplet_attributes_1d{FT}(ξstart, volumes,dry_r3,z_loc, cell_id,grid_range)
    return droplets
end

function plot_env_profiles(grid)
    p1 = plot(grid.states.θ,grid.centers_z, ylabel="Height (m)", xlabel="Potential Temperature", title="θ", legend=false)
    p2 = plot(grid.states.qv*1000,grid.centers_z, ylabel="Height (m)", xlabel="q_tot (g/kg)", title="q_v", legend=false)
    p6 = plot((grid.diagnostics.cloud_LWC+grid.diagnostics.rain_LWC)*1000,grid.centers_z, ylabel="Height (m)", xlabel="LWC", title="LWC", legend=false)
    # p6 = plot((grid.diagnostics.aerosol_LWC)*1000,grid.centers_z, ylabel="Height (m)", xlabel="LWC", title="LWC", legend=false)

    p3 = plot(grid.states.P,grid.centers_z, ylabel="Height (m)", xlabel="Pressure (Pa)", title="P", legend=false)
    p4 = plot(grid.states.ρ,grid.centers_z, ylabel="Height (m)", xlabel="density (kg/m3)", title="ρ", label = "ρ")
    # p4 = vline!([ρ_inv], line=:dash, color=:red, label="ρ_inv")
    p5 = plot(grid.states.T,grid.centers_z, ylabel="Height (m)", xlabel="Temperature(K)", title="T", legend=false)
    p5 = plot!(T_virtual.(grid.states.T, grid.states.qv),grid.centers_z, ylabel="Height (m)", xlabel="Temperature(K)", title="T", legend=false)
    plot(p1, p2, p3,p4,p5,p6, layout=(3,2), size=(800,900))
end


function create_condensation_integrator(grid, drops, condensationsettings, coagsettings, spatialsettings,constants)
    nz = grid.nz
    lnR = log.(volume_to_radius.(drops.X))
    Y = ComponentVector{FT}(lnR = lnR, qvap=grid.states.qv, T = grid.states.T)
    p = (nz,drops, grid.states, constants, condensationsettings,spatialsettings)
    condensation_prob = ODEProblem(condensation_rhs!, Y, (0.0, 1.0), p)
    # condensation_prob = ODEProblem(condensation_rhs_single_cell, Y, (0.0, 1.0), p)

    condensation_integrator = init(condensation_prob, ImplicitEuler(), reltol = 1e-12, abstol = 1e-12)
    return condensation_integrator
end