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


function initialize_scm_environment(nz, dz, P_surface, θl, qt, prescribed_w,init_dist,coagsettings,spatialsettings,condensationsettings,tkesettings,constants)
    droplets = init_droplets_dycoms_scm(init_dist,coagsettings, spatialsettings)
    grid = create_scm_grids(nz, dz,droplets,spatial = spatialsettings,output=nothing)

    p = (constants, θl, qt)

    press = ODEProblem(dP_dz,P_surface,(0,maximum(grid.centers_z)+1),p)
    z_save = sort(unique([grid.centers_z; grid.faces_z]))
    P_init = solve(press,Tsit5(), reltol = 1e-10, abstol = 1e-10,saveat = z_save, tstops = [740.0])
    grid.states.P .= P_init.(grid.centers_z)
    grid.states.P_faces .= P_init.(grid.faces_z)

    grid.states.θ .= θl.(grid.centers_z)
    grid.states.qv .= qt.(grid.centers_z)
    # grid.states.T .= T_from_theta(grid.states.θ,grid.states.P,constants)
    # grid.states.ρ .= grid.states.P ./ (constants.Rd .* T_virtual.(grid.states.T, grid.states.qv))
    # grid.states.ρ .= ρ_ideal_gas(grid.states.P, grid.states.T, grid.states.qv, constants)
    grid.states.ρ .= ρ_calc_θ.(grid.states.P,grid.states.θ,grid.states.qv,constants)
    ρ_dry = grid.states.ρ ./ (1 .+ grid.states.qv)

    # for k in 1:nz
        
    #     grid.states.e[k] = grid.states.qv[k] > 0.008 ? 0.01 : 0
    # end

    grid.wind.u .= tkesettings.geostrophic_u.(grid.centers_z)
    grid.wind.v .= tkesettings.geostrophic_v.(grid.centers_z)
    grid.wind.w .= prescribed_w.(grid.faces_z)

        #set to eq radius, rescale multiplicities from STP to actual density
    FT = eltype(grid.states.P)
    ρ_STP = FT(101325 / (constants.Rd * (273.15)))
    for k in range(1,nz)
        drop_idx = findall(i -> droplets.cell_id[i] == k, 1:length(droplets.X))
        ρ_ratio = ρ_dry[k] / ρ_STP
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
    coagdata = coagulation_run_spatial{FT}(nz, coagsettings.Ns,droplets)
    # condensation_integrator = create_condensation_integrator(grid, droplets, condensationsettings, coagsettings, spatialsettings,constants)
    conddata = condensation_data(FT,nz,coagsettings.Ns)
    raddata = radiation_data(FT,length(droplets.X),grid,constants)
    mpdatatmp = mpdata_tmp_1d(grid.states.qv, grid.faces_z)
    turbdata = turbulence_data(FT, nz)
    return grid,droplets,coagdata,conddata,raddata,mpdatatmp,turbdata
end


function init_droplets_dycoms_scm(dist,settings::coag_settings{FT},
    spatial::spatial_settings{FT})where FT<:AbstractFloat

    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0



    percentile_limit = [0.001, 0.999]  
    Rmin = quantile(dist, percentile_limit[1])
    Rmax = quantile(dist, percentile_limit[2])   
    ξstart::Vector{Int} = fill(1, Ns) 
    z_loc = fill(FT(0.0), Ns) #initialize all droplets at middle of domain
    rad = fill(FT(0.0), Ns) #initialize as just aerosol volume.
    cell_id = fill(1, Ns) #initialize as just aerosol volume.

    # Rarray = range(Rmin, Rmax, length=Ns+1)
    # cdf_values = cdf.(dist, Rarray)
    # rad = (Rarray[2:end] .+ Rarray[1:end-1])./ 2 
    # cdf_values = cdf_values[2:end] - cdf_values[1:end-1]
    nz = spatial.Nz
    Ns_per_grid = Int(floor(Ns / nz))
    for k in 1:nz
        drop_idx = (k-1)*Ns_per_grid+1 : min(k*Ns_per_grid, Ns)
        
        cdf_values = sort(rand(Uniform(),2*Ns_per_grid+1))
        rad[drop_idx] .= quantile.(dist, cdf_values[2:2:end-1])
        cdf_values = cdf_values[1:2:end]
        cdf_vals = cdf_values[2:end] - cdf_values[1:end-1]

        multiplicities = cdf_vals .* (n0 * ΔV)
        ξstart[drop_idx] .= floor.(Int, multiplicities[1:length(drop_idx)] .+ 0.5)
        z_loc[drop_idx] .= (k-1)*spatial.z_grid_height .+ spatial.z_grid_height.*rand(Ns_per_grid)
                
        cell_id[drop_idx] .= k
    end

    dry_r3 = rad.^3
    volumes = 4* pi / 3 .* dry_r3 #initialize as just aerosol volume.

    dropidx = collect(1:Ns)

    sort!(dropidx, by = i -> cell_id[i])
    grid_range = map(1:spatial.Nz) do g
        s = findfirst(i -> cell_id[i] == g, dropidx)
        s === nothing ? (1:0) : s:findlast(i -> cell_id[i] == g, dropidx)
    end
    w_prime = zeros(FT, Ns) 

    droplets = droplet_attributes_1d{FT}(ξstart, volumes,dry_r3,z_loc, cell_id,w_prime,grid_range, dropidx)
    return droplets
end

function plot_env_profiles(grid)
    compute_ql_at_cell!.(grid.states, 1:grid.nz,constants)
    ql = grid.states.ql_tmp
    grid.states.T_tmp .= T_from_theta.(grid.states.θ, grid.states.P, constants)
    p1 = plot(grid.states.θ,grid.centers_z, ylabel="Height (m)", xlabel="Potential Temperature", title="θ", label="θ")
    p1 = plot!(θl.(grid.states.P,grid.states.T_tmp,ql,constants),grid.centers_z, ylabel="Height (m)", xlabel="Potential Temperature", label="θl")
    p2 = plot(grid.states.qv*1000,grid.centers_z, ylabel="Height (m)", xlabel="q (g/kg)", title="q_v", label="qv")
    p2 = plot!(ql*1000,grid.centers_z, ylabel="Height (m)", xlabel="q (g/kg)", title="q_v", legend=false, label="ql")
    p2 = plot!(grid.states.qv*1000 .+ ql*1000,grid.centers_z, ylabel="Height (m)", xlabel="q (g/kg)", title="q_v", legend=false, label="qt")
    p6 = plot(ql.*1000,grid.centers_z, ylabel="Height (m)", xlabel="ql", title="ql",label="tot")
    # p6 = plot!((grid.diagnostics.aerosol_LWC)./grid.states.ρ*1000,grid.centers_z, ylabel="Height (m)", xlabel="ql", title="ql", label="Aerosol")
    # p6 = plot!((grid.diagnostics.cloud_LWC)./grid.states.ρ*1000,grid.centers_z, ylabel="Height (m)", xlabel="ql", title="ql", label="Cloud")
    p8 = plot((grid.diagnostics.rain_LWC)./grid.states.ρ*1000,grid.centers_z, ylabel="Height (m)", xlabel="rain_ql", title="rain_ql", label="Rain")

    p3 = plot(grid.states.P,grid.centers_z, ylabel="Height (m)", xlabel="Pressure (Pa)", title="P", legend=false)
    p4 = plot(grid.states.ρ,grid.centers_z, ylabel="Height (m)", xlabel="density (kg/m3)", title="ρ", label = "ρ")
    # p4 = vline!([ρ_inv], line=:dash, color=:red, label="ρ_inv")
    T_from_theta!(grid.states.T_tmp,grid.states.θ,grid.states.P,constants)
    p5 = plot(grid.states.T_tmp,grid.centers_z, ylabel="Height (m)", xlabel="Temperature(K)", title="T", legend=false)
    p5 = plot!(T_virtual.(grid.states.T_tmp, grid.states.qv),grid.centers_z, ylabel="Height (m)", xlabel="Temperature(K)", title="T", legend=false)
    
    p7 = plot(grid.states.e,grid.centers_z, ylabel="Height (m)", xlabel="Turbulent Kinetic Energy (J/kg)", title="TKE", legend=false)

    envplot = plot(p1, p2, p3,p4,p5,p6,p7,p8, layout=(4,2), size=(800,900))
    return envplot
end

function plot_output_timeseries(grid)
    time = grid.output.time ./ 3600
    t_skips = Int(floor(length(time)/5))
    z_centers = grid.centers_z
    p1 = plot(grid.output.θ[:,1:t_skips:end], z_centers, ylabel="z (m)", title="Potential Temperature", xlabel="θ", legend=false)
    p2 = plot(grid.output.qv[:,1:t_skips:end]*1000, z_centers, ylabel="z (m)", xlabel="qv (g/kg)", title="qv", legend=false)
    p3 = plot(grid.output.cloud_heating_rate[:,1:t_skips:end]*3600, z_centers, ylabel="z (m)", xlabel="Heating Rate (K/hr)", title="Heating Rate", legend=false)
    LWP = grid.output.LWP
    p4 = plot(time, LWP*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries", label = "Total")
    LWP_cloud = sum(grid.output.cloud_LWC, dims=1)' .* grid.dz
    p4 = plot!(time, LWP_cloud*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries", label="Cloud")
    LWP_rain = sum(grid.output.rain_LWC, dims=1)' .* grid.dz
    p4 = plot!(time, LWP_rain*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries", label="Rain")
    # LWP_aerosol = sum(grid.output.aerosol_LWC, dims=1)' .* grid.dz
    # p4 = plot!(time, LWP_aerosol*1000, xlabel="Time (h)", ylabel="LWP (g/m2)", title="LWP Timeseries",label="Aerosol")
    p5 = plot(time, grid.output.surface_precipitation *3600 * 24 * 1000, xlabel="Time (h)", ylabel="pr (mm/day)", title="Surface Precipitation", legend=false)
    # p6 = plot((grid.output.cloud_LWC./grid.output.ρ)[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="Cloud ql (g/kg)", title="Cloud Liquid (not rain)", legend=false)
    p6 = plot(grid.output.ql[:,1:t_skips:end]*1000, z_centers, ylabel="z (m)", xlabel="ql (g/kg)", title="ql", legend=false)
    # p6 = plot!(grid.output.rain_LWC[:,1:5:end]*1000, z_centers, ylabel="z (m)", xlabel="Rain LWC (g/kg)", title="Rain LWC", legend=false)
    
    θl_z = θl.(grid.output.P[:,1:t_skips:end],T_from_theta.(grid.output.θ[:,1:t_skips:end],grid.output.P[:,1:t_skips:end],constants), grid.output.ql[:,1:t_skips:end], constants)
    p7 = plot(θl_z, z_centers, ylabel="z (m)", title="Liquid Water Potential Temperature", xlabel="θl", legend=false)
    p8 = plot(grid.output.ql[:,1:t_skips:end]*1000 .+ grid.output.qv[:,1:t_skips:end]*1000, z_centers, ylabel="z (m)", xlabel="qt(g/kg)", title="qtot", legend=false)
    p9 = plot(grid.output.e[:,1:t_skips:end], z_centers, ylabel="z (m)", xlabel="TKE [m2/s2]", title="TKE", legend=false)

    #p10 is inversion height and cloud base height timeseries. these need to be calculated from the profiles
    #inversion is where qt = 8 g/kg and cloud base is where cloud LWC > 0.01 g/m3
    inv_height = zeros(length(time))
    cloud_base_height = zeros(length(time))
    for i in 1:length(time)
        qt_profile = grid.output.ql[:,i] .+ grid.output.qv[:,i]
        inv_idx = findfirst(qt_profile .< 0.008) 
        inv_height[i] = z_centers[inv_idx]
        # cloud_lwc_profile = grid.output.cloud_LWC[:,i]
        # cloud_base_idx = findfirst(cloud_lwc_profile[5:end] .> 0.1/1000)
        # cloud_base_height[i] = cloud_base_idx == nothing ? NaN : z_centers[cloud_base_idx+4]        
        n_prof = grid.output.number[:,i]
        cloud_base_idx = findfirst(n_prof[5:end] .> 0)
        cloud_base_height[i] = cloud_base_idx == nothing ? NaN : z_centers[cloud_base_idx]
    end
    p10 = plot(time, inv_height, xlabel="Time (h)", ylabel="Inversion Height (m)", label="z_inv")
    p10 = plot!(time, cloud_base_height, xlabel="Time (h)", ylabel="Height (m)", label="z_cb")
    # p11 = plot(grid.output.condensation_rad_net[:,1:t_skips:end]*1000, z_centers, ylabel="Height (m)", xlabel="dqv[g/kg•s]", title="rad cond dqv", legend=false)
    meancloudnumber = [let col = filter(>(0), grid.output.number[:,t]); isempty(col) ? 0.0 : mean(col); end for t in axes(grid.output.number, 2)]
    p11 = plot(time,meancloudnumber*1e-6/50, ylabel="n_c [cm-3]", xlabel="", title="Time (h)", legend=false)
    
    condensation_time_series = sum(grid.output.condensation_src, dims=1)'
    condensation_time_series_rad_net = sum(grid.output.condensation_rad_net, dims=1)'
    condensation_time_series_rad_abs = sum(grid.output.condensation_rad_abs, dims=1)'
    if maximum(time) >1
        spinuptimeidx = findfirst(time .> 1)
        p12 = plot(time[spinuptimeidx:end],condensation_time_series[spinuptimeidx:end]*3600*1000, xlabel="Time (h)", ylabel="dqv[g/kg•hr]", title="cond dqv", label="total")
        p12 = plot!(time[spinuptimeidx:end],condensation_time_series_rad_net[spinuptimeidx:end]*3600*1000, xlabel="Time (h)", ylabel="dqv[g/kg•hr]", title="cond dqv", label="radiation")
        # p12 = plot!(time[spinuptimeidx:end],condensation_time_series_rad_abs[spinuptimeidx:end]*3600*1000, xlabel="Time (h)", ylabel="dqv[g/kg•hr]", title="cond dqv", label="Rad Abs")
    else
        p12 = plot(time,condensation_time_series*3600, xlabel="Time (h)", ylabel="Condensation Mass Source (kg/m3/hr)", title="Condensation Mass Source", label="Mass")
    end
    tplot = plot(p1, p2,p7,p8, p3, p4,p5,p6,p9,p10,p11,p12, layout=(4,3), size=(900,900))
    return tplot
end

function plot_ensemble_timeseries(rad_output, base_output, ref_grid)
    seeds_rad  = collect(keys(rad_output))
    seeds_base = collect(keys(base_output))
    time       = base_output[first(seeds_base)].time ./ 3600
    z_centers  = ref_grid.centers_z
    n_t        = length(time)
    t_skips    = max(1, Int(floor(n_t / 3)))
    t_idx      = 1:t_skips:n_t

    function ts_ribbon(output_dict, seeds, field)
        mat = hcat([getfield(output_dict[s], field) for s in seeds]...)
        med = [median(mat[i, :]) for i in 1:n_t]
        lo  = [quantile(mat[i, :], 0.25) for i in 1:n_t]
        hi  = [quantile(mat[i, :], 0.75) for i in 1:n_t]
        med, med .- lo, hi .- med
    end

    function profile_stats(output_dict, seeds, field)
        arr = cat([getfield(output_dict[s], field)[:, t_idx] for s in seeds]..., dims=3)
        med = [median(arr[z, i, :]) for z in axes(arr, 1), i in axes(arr, 2)]
        lo  = [quantile(arr[z, i, :], 0.25) for z in axes(arr, 1), i in axes(arr, 2)]
        hi  = [quantile(arr[z, i, :], 0.75) for z in axes(arr, 1), i in axes(arr, 2)]
        med, lo, hi
    end

    function add_ribbon_profiles!(p, med, lo, hi; scale=1.0)
        for i in axes(med, 2)
            shade_x = vcat(lo[:, i] .* scale, reverse(hi[:, i] .* scale))
            shade_y = vcat(z_centers, reverse(z_centers))
            plot!(p, shade_x, shade_y; seriestype=:shape, alpha=0.2, linewidth=0, color=i)
            plot!(p, med[:, i] .* scale, z_centers; color=i, label="")
        end
    end

    function add_dashed_profiles!(p, med; scale=1.0)
        for i in axes(med, 2)
            plot!(p, med[:, i] .* scale, z_centers; color=i, linestyle=:dash, label="")
        end
    end

    function make_profile_plot(field, xlabel, title; scale=1.0)
        med_b, lo_b, hi_b = profile_stats(base_output, seeds_base, field)
        med_r = profile_stats(rad_output, seeds_rad, field)[1]
        p = plot(; title, xlabel, ylabel="z (m)", legend=false)
        add_ribbon_profiles!(p, med_b, lo_b, hi_b; scale)
        add_dashed_profiles!(p, med_r; scale)
        p
    end

    p1 = make_profile_plot(:θ,  "θ (K)",       "Potential Temperature")
    p2 = make_profile_plot(:qv, "qv (g/kg)",   "qv"; scale=1000.0)
    p6 = make_profile_plot(:ql, "ql (g/kg)",   "ql"; scale=1000.0)
    p9 = make_profile_plot(:e,  "TKE [m²/s²]", "TKE")

    function θl_profile_stats(output_dict, seeds)
        arr = cat([θl.(output_dict[s].P[:, t_idx],
                       T_from_theta.(output_dict[s].θ[:, t_idx],
                                     output_dict[s].P[:, t_idx], constants),
                       output_dict[s].ql[:, t_idx], constants)
                   for s in seeds]..., dims=3)
        ([median(arr[z, i, :])         for z in axes(arr,1), i in axes(arr,2)],
         [quantile(arr[z, i, :], 0.25) for z in axes(arr,1), i in axes(arr,2)],
         [quantile(arr[z, i, :], 0.75) for z in axes(arr,1), i in axes(arr,2)])
    end
    θl_b = θl_profile_stats(base_output, seeds_base)
    θl_r = θl_profile_stats(rad_output, seeds_rad)
    p7 = plot(; title="θl", xlabel="θl (K)", ylabel="z (m)", legend=false)
    add_ribbon_profiles!(p7, θl_b[1], θl_b[2], θl_b[3])
    add_dashed_profiles!(p7, θl_r[1])

    function ts_plot(field; scale=1.0, ylabel="", title="")
        med_b, lo_b, hi_b = ts_ribbon(base_output, seeds_base, field)
        med_r, lo_r, hi_r = ts_ribbon(rad_output,  seeds_rad,  field)
        p = plot(time, med_b .* scale; ribbon=(lo_b .* scale, hi_b .* scale),
                 fillalpha=0.3, label="base", ylabel, title, xlabel="Time (h)")
        plot!(p, time, med_r .* scale; ribbon=(lo_r .* scale, hi_r .* scale),
              fillalpha=0.3, label="rad")
        p
    end
    p4 = ts_plot(:LWP;                  scale=1000.0,    ylabel="LWP (g/m²)",   title="LWP")
    p5 = ts_plot(:surface_precipitation; scale=3600*24.0, ylabel="pr (mm/day)", title="Precipitation")

    function cloud_layer_heights(output_dict, seeds)
        qt_arr  = cat([output_dict[s].ql  .+ output_dict[s].qv for s in seeds]..., dims=3)
        num_arr = cat([output_dict[s].number for s in seeds]..., dims=3)
        inv_h = zeros(n_t); cb_h = fill(NaN, n_t)
        for t in 1:n_t
            qt_med  = [median(qt_arr[z,  t, :]) for z in axes(qt_arr,  1)]
            num_med = [median(num_arr[z, t, :]) for z in axes(num_arr, 1)]
            idx_inv = findfirst(qt_med .< 0.008)
            inv_h[t] = idx_inv === nothing ? NaN : z_centers[idx_inv]
            idx_cb  = findfirst(num_med[5:end] .> 0)
            cb_h[t]  = idx_cb  === nothing ? NaN : z_centers[idx_cb]
        end
        inv_h, cb_h
    end
    inv_b, cb_b = cloud_layer_heights(base_output, seeds_base)
    inv_r, cb_r = cloud_layer_heights(rad_output,  seeds_rad)

    p10 = plot(time, inv_b; label="z_inv",      ylabel="Height (m)", xlabel="Time (h)", title="Cloud Layer Heights")
    plot!(p10, time, cb_b;  label="z_cb",       linestyle=:dash)
    plot!(p10, time, inv_r; label="z_inv (rad)", linestyle=:dot)
    plot!(p10, time, cb_r;  label="z_cb (rad)",  linestyle=:dashdot)

    plot(p1, p2, p7, p6, p4, p5, p9, p10; layout=(2, 4), size=(1200, 600))
end

function plot_ensemble_field(field, ensemble_output, ref_grid=nothing;
                            scale=1.0, show_ribbon=true, axis=:height, t_window=nothing, color=:auto, label="", key_filter=nothing, kwargs...)
    all_keys = collect(keys(ensemble_output))
    # filter by non-seed key prefix when keys are tuples
    selected = key_filter === nothing ? all_keys :
               filter(k -> k isa Tuple && k[1:end-1] == key_filter, all_keys)
    # drop failed runs
    selected = filter(k -> ensemble_output[k] !== nothing, selected)

    sample = getfield(ensemble_output[first(selected)], field)
    time   = ensemble_output[first(selected)].time ./ 3600

    if ndims(sample) == 1
        mat = hcat([getfield(ensemble_output[s], field) for s in selected]...)
        mn  = vec(mean(mat, dims=2))
        sd  = vec(std(mat,  dims=2))
        rib = show_ribbon ? sd .* scale : nothing
        plot!(time, mn .* scale; ribbon=rib, fillalpha=0.3, xlabel="Time (h)", label=label, color=color, kwargs...)
    elseif axis == :time
        mat = hcat([vec(mean(getfield(ensemble_output[s], field), dims=1)) for s in selected]...)
        mn  = vec(mean(mat, dims=2))
        sd  = vec(std(mat,  dims=2))
        rib = show_ribbon ? sd .* scale : nothing
        plot!(time, mn .* scale; ribbon=rib, fillalpha=0.3, xlabel="Time (h)", label=label, color=color, kwargs...)
    else  # axis == :height
        z_centers = ref_grid.centers_z
        t_idx = t_window === nothing ? axes(sample, 2) :
                findall(t -> t_window[1] <= time[t] <= t_window[2], axes(sample, 2))
        mat = hcat([vec(mean(getfield(ensemble_output[s], field)[:, t_idx], dims=2)) for s in selected]...)
        mn  = vec(mean(mat, dims=2))
        sd  = vec(std(mat,  dims=2))
        p = plot!(; ylabel="z (m)", kwargs...)
        if show_ribbon
            shade_x = vcat((mn .- sd) .* scale, reverse((mn .+ sd) .* scale))
            shade_y = vcat(z_centers, reverse(z_centers))
            plot!(p, shade_x, shade_y; seriestype=:shape, alpha=0.2, linewidth=0, color=color, label="")
        end
        plot!(p, mn .* scale, z_centers; color=color, label=label, kwargs...)
        p
    end
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