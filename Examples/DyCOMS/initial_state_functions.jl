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


function initialize_scm_environment(nz, dz, P_surface, θl, qt, geostrophic_u, geostrophic_v, prescribed_w,droplets,spatialsettings,condensationsettings,constants)
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

        #set to eq radius, rescale multiplicities from STP to actual density
    FT = eltype(grid.states.P)
    ρ_STP = FT(101325 / (constants.Rd * 273.15))
    for k in range(1,nz)
        drop_idx = findall(i -> droplets.cell_id[i] == k, 1:length(droplets.X))
        ρ_ratio = grid.states.ρ[k] / ρ_STP
        droplets.ξ[drop_idx] .= floor.(Int, droplets.ξ[drop_idx] .* ρ_ratio .+ 0.5)

        T = T_from_theta(grid.states.θ[k], grid.states.P[k], constants)
        qv_k = grid.states.qv[k]
        P_k = grid.states.P[k]
        S_env = sat(qv_k,P_k)/esat(T)
        if S_env > 0.95
            S_env = 0.95
        end
        # println(S_env)
        find_equilibrium_radius.(droplets,drop_idx, condensationsettings.kappa, T, S_env,constants)
        
    end
    sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)
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
    w_prime = zeros(FT, Ns) 

    droplets = droplet_attributes_1d{FT}(ξstart, volumes,dry_r3,z_loc, cell_id,w_prime,grid_range, dropidx)
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

    envplot = plot(p1, p2, p3,p4,p5,p6,p7, layout=(4,2), size=(800,900))
    return envplot
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
    p5 = plot(time, grid.output.surface_precipitation *3600 * 24, xlabel="Time (h)", ylabel="pr (mm/day)", title="Surface Precipitation", legend=false)
    # p6 = plot(grid.output.cloud_LWC[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="Cloud LWC (g/m3)", title="Cloud LWC", legend=false)
    p6 = plot(grid.output.ql[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="ql (g/kg)", title="Cloud LWC", legend=false)
    # p6 = plot!(grid.output.rain_LWC[:,1:5:end]*1000, z_centers, ylabel="Height (m)", xlabel="Rain LWC (g/kg)", title="Rain LWC", legend=false)
    
    θl = θl_θ.(grid.output.θ[:,1:t_skips:end], grid.output.ql[:,1:t_skips:end], constants)
    p7 = plot(θl, z_centers, ylabel="Height (m)", xlabel="Liquid Water Potential Temperature", title="θl", legend=false)
    p8 = plot(grid.output.ql[:,1:t_skips:end]*1000 .+ grid.output.qv[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="q tot(g/kg)", title="qtot", legend=false)
    p9 = plot(grid.output.e[:,1:t_skips:end], z_centers, ylabel="Height (m)", xlabel="Turbulent Kinetic Energy (J/kg)", title="TKE", legend=false)

    #p10 is inversion height and cloud base height timeseries. these need to be calculated from the profiles
    #inversion is where qt = 8 g/kg and cloud base is where cloud LWC > 0.01 g/m3
    inv_height = zeros(length(time))
    cloud_base_height = zeros(length(time))
    for i in 1:length(time)
        qt_profile = grid.output.ql[:,i] .+ grid.output.qv[:,i]
        inv_idx = findfirst(qt_profile .< 0.008)
        inv_height[i] = z_centers[inv_idx]
        cloud_lwc_profile = grid.output.cloud_LWC[:,i]
        cloud_base_idx = findfirst(cloud_lwc_profile[5:end] .> 0.1/1000)
        cloud_base_height[i] = cloud_base_idx == nothing ? NaN : z_centers[cloud_base_idx+4]
    end
    p10 = plot(time, inv_height, xlabel="Time (h)", ylabel="Inversion Height (m)", label="z_inv")
    p10 = plot!(time, cloud_base_height, xlabel="Time (h)", ylabel="Height (m)", label="z_cb")
    p11 = plot(grid.output.condensation_rad_abs[:,1:t_skips:end], z_centers, ylabel="Height (m)", xlabel="Condensation Rad Abs (W/m3)", title="Condensation Rad Abs", legend=false)
    condensation_time_series = sum(grid.output.condensation_src, dims=1)'
    condensation_time_series_rad_net = sum(grid.output.condensation_rad_net, dims=1)'
    condensation_time_series_rad_abs = sum(grid.output.condensation_rad_abs, dims=1)'
    if maximum(time) >1
        spinuptimeidx = findfirst(time .> 1)
        p12 = plot(time[spinuptimeidx:end],condensation_time_series[spinuptimeidx:end]*3600, xlabel="Time (h)", ylabel="Condensation Mass Source (kg/m3/hr)", title="Condensation Mass Source", label="Mass")
        p12 = plot!(time[spinuptimeidx:end],condensation_time_series_rad_net[spinuptimeidx:end]*3600, xlabel="Time (h)", ylabel="Condensation Rad Net (W/m3)", title="Condensation Mass Source", label="Rad Net")
        p12 = plot!(time[spinuptimeidx:end],condensation_time_series_rad_abs[spinuptimeidx:end]*3600, xlabel="Time (h)", ylabel="Condensation Rad Abs (W/m3)", title="Condensation Mass Source", label="Rad Abs")
    else
        p12 = plot(time,condensation_time_series*3600, xlabel="Time (h)", ylabel="Condensation Mass Source (kg/m3/hr)", title="Condensation Mass Source", label="Mass")
    end
    tplot = plot(p1, p2,p7,p8, p3, p4,p5,p6,p9,p10,p11,p12, layout=(4,3), size=(900,900))
    return tplot
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