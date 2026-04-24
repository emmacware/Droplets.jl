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


function initialize_scm_environment(nz, dz, P_surface, θl, qt, geostrophic_u, geostrophic_v, prescribed_w,droplets,spatialsettings)
    grid = create_scm_grids(nz, dz,droplets,spatial = spatialsettings,output=nothing)

    p = (constants, θl, qt)

    press = ODEProblem(dP_dz,P_surface,(0,maximum(grid.centers_z)+1),p)
    z_save = sort(unique([grid.centers_z; grid.faces_z]))
    P_init = solve(press,Tsit5(), reltol = 1e-10, abstol = 1e-10,saveat = z_save)
    grid.states.P .= P_init.(grid.centers_z)
    grid.states.P_faces .= P_init.(grid.faces_z)

    grid.states.θ .= θl.(grid.centers_z)
    grid.states.qv .= qt.(grid.centers_z)
    # grid.states.T .= T_from_theta(grid.states.θ,grid.states.P,constants)
    # grid.states.ρ .= grid.states.P ./ (constants.Rd .* T_virtual.(grid.states.T, grid.states.qv))
    # grid.states.ρ .= ρ_ideal_gas(grid.states.P, grid.states.T, grid.states.qv, constants)
    grid.states.ρ .= ρ_calc_θ(grid.states.P,grid.states.θ,grid.states.qv,constants)

    grid.wind.u .= geostrophic_u.(grid.centers_z)
    grid.wind.v .= geostrophic_v.(grid.centers_z)
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


    droplets = droplet_attributes_1d{FT}(ξstart, volumes,dry_r3,z_loc, cell_id,grid_range, dropidx)
    return droplets
end

function plot_env_profiles(grid)
    ql = compute_ql_at_cell.(grid.states, collect(1:grid.nz))
    p1 = plot(grid.states.θ,grid.centers_z, ylabel="Height (m)", xlabel="Potential Temperature", title="θ", label="θ")
    p1 = plot!(θl_θ.(grid.states.θ,ql,constants),grid.centers_z, ylabel="Height (m)", xlabel="Potential Temperature", label="θl")
    p2 = plot(grid.states.qv*1000,grid.centers_z, ylabel="Height (m)", xlabel="q (g/kg)", title="q_v", label="qv")
    p2 = plot!(ql*1000,grid.centers_z, ylabel="Height (m)", xlabel="q (g/kg)", title="q_v", legend=false, label="ql")
    p2 = plot!(grid.states.qv*1000 .+ ql*1000,grid.centers_z, ylabel="Height (m)", xlabel="q (g/kg)", title="q_v", legend=false, label="qt")
    p6 = plot((grid.diagnostics.cloud_LWC+grid.diagnostics.rain_LWC+grid.diagnostics.aerosol_LWC)*1000,grid.centers_z, ylabel="Height (m)", xlabel="LWC", title="LWC")
    p6 = plot!((grid.diagnostics.aerosol_LWC)*1000,grid.centers_z, ylabel="Height (m)", xlabel="LWC", title="LWC", label="Aerosol")
    p6 = plot!((grid.diagnostics.cloud_LWC)*1000,grid.centers_z, ylabel="Height (m)", xlabel="LWC", title="LWC", label="Cloud")
    p6 = plot!((grid.diagnostics.rain_LWC)*1000,grid.centers_z, ylabel="Height (m)", xlabel="LWC", title="LWC", label="Rain")

    p3 = plot(grid.states.P,grid.centers_z, ylabel="Height (m)", xlabel="Pressure (Pa)", title="P", legend=false)
    p4 = plot(grid.states.ρ,grid.centers_z, ylabel="Height (m)", xlabel="density (kg/m3)", title="ρ", label = "ρ")
    # p4 = vline!([ρ_inv], line=:dash, color=:red, label="ρ_inv")
    grid.states.T_tmp .= T_from_theta.(grid.states.θ,grid.states.P,constants)
    p5 = plot(grid.states.T_tmp,grid.centers_z, ylabel="Height (m)", xlabel="Temperature(K)", title="T", legend=false)
    p5 = plot!(T_virtual.(grid.states.T_tmp, grid.states.qv),grid.centers_z, ylabel="Height (m)", xlabel="Temperature(K)", title="T", legend=false)
    
    p7 = plot(grid.states.e,grid.centers_z, ylabel="Height (m)", xlabel="Turbulent Kinetic Energy (J/kg)", title="TKE", legend=false)

    plot(p1, p2, p3,p4,p5,p6,p7, layout=(4,2), size=(800,900))
end

function plot_output_timeseries(grid)
    time = grid.output.time ./ 3600
    t_skips = Int(floor(length(time)/5))
    z_centers = grid.centers_z
    p1 = plot(grid.output.θ[:,1:t_skips:end], z_centers, ylabel="Height (m)", xlabel="Potential Temperature", title="θ", legend=false)
    p2 = plot(grid.output.qv[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="qv (g/kg)", title="qv", legend=false)
    p3 = plot(grid.output.cloud_heating_rate[:,1:t_skips:end]*3600, z_centers, ylabel="Height (m)", xlabel="Cloud Heating Rate (K/hr)", title="Cloud Heating Rate", legend=false)
    LWP = grid.output.LWP
    p4 = plot(time, LWP*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries", label = "Total")
    LWP_cloud = sum(grid.output.cloud_LWC, dims=1)' .* grid.dz
    p4 = plot!(time, LWP_cloud*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries", label="Cloud")
    LWP_rain = sum(grid.output.rain_LWC, dims=1)' .* grid.dz
    p4 = plot!(time, LWP_rain*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries", label="Rain")
    LWP_aerosol = sum(grid.output.aerosol_LWC, dims=1)' .* grid.dz
    p4 = plot!(time, LWP_aerosol*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries",label="Aerosol")
    p5 = plot(time, grid.output.surface_precipitation *3600, xlabel="Time (h)", ylabel="Surface Precipitation (mm/hr)", title="Surface Precipitation Timeseries", legend=false)
    p6 = plot(grid.output.cloud_LWC[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="Cloud LWC (g/m3)", title="Cloud LWC", legend=false)
    # p6 = plot!(grid.output.rain_LWC[:,1:5:end]*1000, z_centers, ylabel="Height (m)", xlabel="Rain LWC (g/kg)", title="Rain LWC", legend=false)
    
    θl = θl_θ.(grid.output.θ[:,1:t_skips:end], grid.output.ql[:,1:t_skips:end], constants)
    p7 = plot(θl, z_centers, ylabel="Height (m)", xlabel="Liquid Water Potential Temperature", title="θl", legend=false)
    p8 = plot(grid.output.ql[:,1:t_skips:end]*1000 .+ grid.output.qv[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="q tot(g/kg)", title="qtot", legend=false)
    plot(p1, p2,p7,p8, p3, p4,p5,p6, layout=(4,2), size=(600,900))
end
# function create_condensation_integrator(grid, drops, condensationsettings, coagsettings, spatialsettings,constants)
#     nz = grid.nz
#     lnr = log.(volume_to_radius.(drops.X))
#     Y = [(lnr[1])]
#     p = (1,drops, grid.states, constants, condensationsettings,spatialsettings)
#     condensation_prob = ODEProblem(condensation_rhs_single_droplet, Y, (0.0, 0.01), p)
#     # condensation_prob = ODEProblem(condensation_rhs_single_cell, Y, (0.0, 1.0), p)

#     condensation_integrator = init(condensation_prob, Rosenbrock23(), reltol = 1e-12, abstol = 1e-12)
#     return condensation_integrator
# end


# p = (1,droplets, grid.states, constants, condensationsettings,spatialsettings)

# du = zeros(1)
# u = [((volume_to_radius.(droplets.X))[1])^2]
# condensation_rhs_single_droplet(du,u,p,0.0)