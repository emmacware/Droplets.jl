using RRTMGP
using NCDatasets

using RRTMGP.LookUpTables

#these inputs have the same names as the scm.jl file so if you run that first you can
#call this function without it needing to be part of the timestepper

function stand_in_radiation_update_function!(grid,spatialsettings, diagnosticsettings, constants)
    
    FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64

    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)

    #number of grid points and grid spacing
    Nz = grid.nz
    dz = grid.dz
    #Face (Nz+1) and center (Nz) heights
    z_faces = grid.faces_z
    z_centers = grid.centers_z

    #Profiles
    θ, qv, T, P, ρ = grid.states.θ, grid.states.qv, grid.states.T, grid.states.P, grid.states.ρ
    # plot_env_profiles(grid) will plot these

    #Liquid water diagnostics (cloud_, aerosol_, and rain_)
    cloud_LWC = grid.diagnostics.cloud_LWC
    cloud_effective_radius = grid.diagnostics.cloud_effective_radius

    lw_file = "/home/aigel/Droplets.jl/data/rrtmgp-data-lw-g256-2018-12-04.nc"
    sw_file = "/home/aigel/Droplets.jl/data/rrtmgp-data-sw-g224-2018-12-04.nc"
    lw_cld_file = "/home/aigel/Droplets.jl/data/rrtmgp-cloud-optics-coeffs-lw.nc"
    sw_cld_file = "/home/aigel/Droplets.jl/data/rrtmgp-cloud-optics-coeffs-reordered-sw.nc"

    # reading longwave gas optics lookup data
    lookup_lw, idx_gases = Dataset(lw_file,"r") do ds 
        LookUpLW(ds, FT, DA)
    end
    
     # reading longwave cloud lookup data ###missing diamice_lwr
     lookup_lw_cld = Dataset(lw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end
    #reading shortwave gas optics lookup data
    lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
        LookUpSW(ds, FT, DA)
    end
    # reading shortwave cloud lookup data ###missing diamice_lwr
    lookup_sw_cld = Dataset(sw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end


    ########This part could be a separate function######
    deg2rad = FT(π) / FT(180)
    nlay = Nz
    #ncol = Int(ds_in.dim["col"]) # col#1 repeated 128 times, per RRTMGP example
    nlev = nlay + 1
    ngas = LookUpTables.get_n_gases(lookup_lw)
    nbnd_lw = LookUpTables.get_n_bnd(lookup_lw)
    nbnd_sw = LookUpTables.get_n_bnd(lookup_sw)
    #---no lat / long information for this
    lon = nothing
    lat = nothing
    # The example only reads the first column and
    # replicates it ncol times.

    sfc_emis = DA{FT, 1}(undef, nbnd_lw)
    sfc_alb_direct = DA{FT, 1}(undef, nbnd_sw)
    sfc_alb_diffuse = DA{FT, 1}(undef, nbnd_sw)
    cos_zenith = DA{FT, 0}(undef)
    irrad = DA{FT, 0}(undef)
    # these values are taken from the example
    sfc_emis .= FT(0.98)
    sfc_alb_direct .= FT(0.06)
    sfc_alb_diffuse .= FT(0.06)
    cos_zenith .= FT(0.86)
    irrad .= FT(lookup_sw.solar_src_tot)


    p_lev = Array{FT}(reshape(Array(ds_in["p_lev"])[1, :], nlev, 1))

    bot_at_1 = p_lev[1, 1] > p_lev[end, 1]
    lev_ind = bot_at_1 ? (1:nlev) : (nlev:-1:1)
    lay_ind = bot_at_1 ? (1:nlay) : (nlay:-1:1)

    p_lev = Array{FT}(reshape(Array(ds_in["p_lev"])[1, lev_ind], nlev, 1))
    p_lay = Array{FT}(reshape(Array(ds_in["p_lay"])[1, lay_ind], nlay, 1))
    t_lev = Array{FT}(reshape(Array(ds_in["t_lev"])[1, lev_ind], nlev, 1))
    t_lay = Array{FT}(reshape(Array(ds_in["t_lay"])[1, lay_ind], nlay, 1))
    t_sfc = Array{FT}(reshape([t_lev[1, 1]], 1))

    #col_dry = DA{FT,2}(transpose(Array(ds_in["col_dry"])[:, lay_ind]))
    #col_dry from the dataset not used in the FORTRAN RRTMGP example

    # Reading volume mixing ratios
    vmrat = zeros(FT, ngas, nlay)

    vmrat[idx_gases["h2o"], :, 1] .= Array{FT}(Array(ds_in["h2o"])[1, lay_ind])
    vmrat[idx_gases["o3"], :, 1] .= Array{FT}(Array(ds_in["o3"])[1, lay_ind])
    vmrat[idx_gases["co2"], :, 1] .= FT(348e-6)
    vmrat[idx_gases["ch4"], :, 1] .= FT(1650e-9)
    vmrat[idx_gases["n2o"], :, 1] .= FT(306e-9)
    vmrat[idx_gases["n2"], :, 1] .= FT(0.7808)
    vmrat[idx_gases["o2"], :, 1] .= FT(0.2095)
    vmrat[idx_gases["co"], :, 1] .= FT(0)

    for icol in 2:ncol
        vmrat[:, :, icol] .= vmrat[:, :, 1]
    end
    vmr = Vmr(DA(vmrat))
    col_dry = DA{FT, 1}(undef, nlay)
    rel_hum = DA{FT, 1}(undef, nlay)
    vmr_h2o = view(vmr.vmr, idx_gases["h2o"], :, :)

    cld_frac = zeros(FT, nlay)
    cld_mask_lw = zeros(Bool, nlay)
    cld_mask_sw = zeros(Bool, nlay)
    cld_r_eff_liq = zeros(FT, nlay)
    cld_r_eff_ice = zeros(FT, nlay)
    cld_path_liq = zeros(FT, nlay)
    cld_path_ice = zeros(FT, nlay)

    radliq_lwr, radliq_upr, radice_lwr, radice_upr = Array(lookup_lw_cld.bounds)
    r_eff_liq = (radliq_lwr + radliq_upr) / FT(2)
    r_eff_ice = (radice_lwr + radice_upr) / FT(2)

    # Restrict clouds to troposphere (> 100 hPa = 100*100 Pa)
    # and not very close to the ground (< 900 hPa), and
    # put them in 2/3 of the columns since that's roughly the
    # total cloudiness of earth
    
    cld_path_liq[ilay, icol] = FT(10)
    cld_r_eff_liq[ilay, icol] = r_eff_liq

    cld_path_ice[ilay, icol] = FT(10)
    cld_r_eff_ice[ilay, icol] = r_eff_ice
 

    p_lay = DA(p_lay)
    p_lev = DA(p_lev)
    t_lay = DA(t_lay)
    t_lev = DA(t_lev)

    device = ClimaComms.device(context)
    compute_col_gas!(device, p_lev, col_dry, param_set, vmr_h2o, lat) # the example skips lat based gravity calculation
    compute_relative_humidity!(device, rel_hum, p_lay, t_lay, param_set, vmr_h2o) # compute relative humidity

    layerdata = similar(p_lay, 4, nlay)
    layerdata[1, :] .= col_dry
    layerdata[2, :] .= p_lay
    layerdata[3, :] .= t_lay
    layerdata[4, :] .= rel_hum

    t_sfc = DA(t_sfc)

    cld_frac = DA(cld_frac)
    cld_mask_lw = DA(cld_mask_lw)
    cld_mask_sw = DA(cld_mask_sw)
    cld_r_eff_liq = DA(cld_r_eff_liq)
    cld_r_eff_ice = DA(cld_r_eff_ice)
    cld_path_liq = DA(cld_path_liq)
    cld_path_ice = DA(cld_path_ice)
    ice_rgh = 2 # medium ice roughness

    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)

    inc_flux = nothing
    slv_lw = SLVLW(grid_params; params = param_set, sfc_emis, inc_flux)
    # Setting up shortwave problem---------------------------------------
    inc_flux_diffuse = nothing
    swbcs = (; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = SLVSW(grid_params; swbcs...)
    #------calling solvers
    metric_scaling = DA(one.(slv_sw.flux.flux_up))
    solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, nothing, metric_scaling)
    if device isa ClimaComms.CPUSingleThreaded
        JET.@test_opt solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, nothing, metric_scaling)
        #@test (@allocated solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, nothing, metric_scaling)) == 0
        @test (@allocated solve_lw!(slv_lw, as, lookup_lw, lookup_lw_cld, nothing, metric_scaling)) ≤ 448
    end

    
    return
end