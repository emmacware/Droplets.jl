using RRTMGP
using NCDatasets
using ClimaComms
using RRTMGP: RRTMGPGridParams
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.RTESolver
using RRTMGP.GrayUtils
import RRTMGP.Parameters.RRTMGPParameters

#these inputs have the same names as the scm.jl file so if you run that first you can
#call this function without it needing to be part of the timestepper

function read_radiation_tables()

    FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)

    lw_file = joinpath(@__DIR__, "../../data/rrtmgp-data-lw-g256-2018-12-04.nc")
    sw_file = joinpath(@__DIR__, "../../data/rrtmgp-data-sw-g224-2018-12-04.nc")
    lw_cld_file = joinpath(@__DIR__, "../../data/rrtmgp-clouds-lw-g256.nc")
    sw_cld_file = joinpath(@__DIR__, "../../data/rrtmgp-clouds-sw-g224.nc")

    # reading longwave gas optics lookup data
    lookup_lw, idx_gases = Dataset(lw_file,"r") do ds 
        LookUpLW(ds, FT, DA)
    end
    
     # reading longwave cloud lookup data 
    lookup_lw_cld = Dataset(lw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end
    #reading shortwave gas optics lookup data
    lookup_sw, idx_gases = Dataset(sw_file, "r") do ds
        LookUpSW(ds, FT, DA)
    end
    # reading shortwave cloud lookup data 
    lookup_sw_cld = Dataset(sw_cld_file, "r") do ds
        LookUpCld(ds, FT, DA)
    end

    return lookup_lw,lookup_lw_cld,lookup_sw,lookup_sw_cld,idx_gases
end

function radiation_function!(grid,spatialsettings, diagnosticsettings, constants, dt, i)
    
    #Some things to consider:
    #1. Do we want incoming downwelling LW?
    #2. Or, do we want to specify the atmosphere above the model top just for radiation?
    #
    #To do:
    #1. Couple heating rates to grid temperature, update other thermo accordingly

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
    

    #overrides = (; grav = 9.80665, molmass_dryair = 0.028964, molmass_water = 0.018016, 
    #            gas_constant = 287.0, kappa_d = 0.28571428571, Stefan = 5.67e-8, avogad = 6.02214076e+23)
    #param_set = RRTMGPParameters(FT,overrides)
    param_set = RRTMGPParameters(9.80665, 0.028964, 0.018016, 287.0, 0.28571428571, 5.67e-8, 6.02214076e+23 )

    ########This part could be a separate function######
    deg2rad = FT(π) / FT(180)
    nlay = Nz
    ncol = 1
    #ncol = Int(ds_in.dim["col"]) # col#1 repeated 128 times, per RRTMGP example
    nlev = nlay + 1
    ngas = LookUpTables.get_n_gases(lookup_lw)
    nbnd_lw = LookUpTables.get_n_bnd(lookup_lw)
    nbnd_sw = LookUpTables.get_n_bnd(lookup_sw)
    #---no lat / long information for this
    lon = nothing
    lat = nothing

    sfc_emis = DA{FT, 1}(undef, nbnd_lw)
    sfc_alb_direct = DA{FT, 1}(undef, nbnd_sw)
    sfc_alb_diffuse = DA{FT, 1}(undef, nbnd_sw)
    cos_zenith = DA{FT, 1}(undef,ncol)
    toa_flux = DA{FT, 1}(undef,ncol)
    # these values are taken from the example
    sfc_emis .= FT(0.98)
    sfc_alb_direct .= FT(0.06)
    sfc_alb_diffuse .= FT(0.06)
    cos_zenith .= FT(0.86)
    toa_flux .= FT(lookup_sw.solar_src_tot)


    #p_lev = Array{FT}(reshape(Array(ds_in["p_lev"])[1, :], nlev, 1))
    p_lay = Array{FT}(reshape(P,nlay,1))

    bot_at_1 = p_lay[1, 1] > p_lay[end, 1]
    lev_ind = bot_at_1 ? (1:nlev) : (nlev:-1:1)
    lay_ind = bot_at_1 ? (1:nlay) : (nlay:-1:1)

    p_lay = Array{FT}(p_lay[lay_ind])
    p_lev = Array{FT}([P_surface;(p_lay[1:nlay-1]+p_lay[2:nlay])/2])
    p_lev = Array{FT}([p_lev;2*p_lay[end]-p_lev[end]]) #could do logarithmic extrapolation instead
    t_lay = Array{FT}(reshape(T[lay_ind],nlay,1))
    t_lev = Array{FT}(t_lay[1:end-1]+t_lay[2:end])/2
    t_lev = Array{FT}([t_lay[1]*2-t_lev[1];t_lev;t_lay[end]*2-t_lev[end]])

    p_lay = Array{FT}(DA(p_lay))
    p_lev = Array{FT}(DA(p_lev))
    t_lay = Array{FT}(DA(t_lay))
    t_lev = Array{FT}(DA(t_lev))

    p_lay = reshape(p_lay,:,1)
    p_lev = reshape(p_lev,:,1)
    t_lay = reshape(t_lay,:,1)
    t_lev = reshape(t_lev,:,1)
   
    t_sfc = Array{FT}(reshape([t_lev[1, 1]], 1))

    # Reading volume mixing ratios
    vmrat = zeros(FT, ngas, nlay, 1)

    vmrat[idx_gases["h2o"], :, 1] .= Array{FT}(reshape(qv[lay_ind] ,nlay,1))
    vmrat[idx_gases["o3"], :, 1] .= FT(0) #no ozone    
    vmrat[idx_gases["co2"], :, 1] .= FT(348e-6)
    vmrat[idx_gases["ch4"], :, 1] .= FT(1650e-9)
    vmrat[idx_gases["n2o"], :, 1] .= FT(306e-9)
    vmrat[idx_gases["n2"], :, 1] .= FT(0.7808)
    vmrat[idx_gases["o2"], :, 1] .= FT(0.2095)
    vmrat[idx_gases["co"], :, 1] .= FT(0)

    vmr = Vmr(DA(vmrat))
    col_dry = DA{FT, 2}(undef, nlay, 1)
    rel_hum = DA{FT, 2}(undef, nlay, 1)
    #vmr_h2o = view(vmr.vmr, idx_gases["h2o"], :, :)
    vmr_h2o = reshape(vmrat[idx_gases["h2o"],:,1],:,1)

    cld_frac = zeros(FT, nlay, 1)
    cld_mask_lw = zeros(Bool, nlay, 1)
    cld_mask_sw = zeros(Bool, nlay, 1)
    cld_r_eff_liq = zeros(FT, nlay, 1)
    cld_r_eff_ice = zeros(FT, nlay, 1)
    cld_path_liq = zeros(FT, nlay, 1)
    cld_path_ice = zeros(FT, nlay, 1)
   
    #Liquid water diagnostics (cloud_, aerosol_, and rain_)
    cld_path_liq .= grid.diagnostics.cloud_LWC * grid.dz * 1000 #units in g/m2
    cld_r_eff_liq .= grid.diagnostics.cloud_effective_radius * 1.0e6 #units in microns
    cld_frac[cld_path_liq.>0.0] .= 1.0
    cld_mask_lw[cld_frac.==1] .= true
    cld_mask_sw[cld_frac.==1] .= true
    
    radliq_lwr = lookup_lw_cld.bounds[1]
    radliq_upr = lookup_lw_cld.bounds[2]

    cld_r_eff_liq[cld_r_eff_liq .< radliq_lwr] .= radliq_lwr
    cld_r_eff_liq[cld_r_eff_liq .> radliq_upr] .= radliq_upr

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
    cloud_state = CloudState(cld_r_eff_liq,cld_r_eff_ice,cld_path_liq,cld_path_ice,cld_frac,cld_mask_lw,cld_mask_sw,MaxRandomOverlap(),ice_rgh)

    grid_params = RRTMGPGridParams(FT; context, nlay, ncol)

    #inc_flux = nothing
    inc_flux = [21.0, 31.0, 5.0, 15.0, 2.0, 0.03, 0.01, 0.04, 1.8, 3.1, 9.3, 0.6, 0.01, 0.4, 0.001, 0.002]
    slv_lw = TwoStreamLWRTE(grid_params; params = param_set, sfc_emis, inc_flux)
    # Setting up shortwave problem---------------------------------------
    inc_flux_diffuse = nothing

    swbcs = (; cos_zenith, toa_flux, sfc_alb_direct, inc_flux_diffuse, sfc_alb_diffuse)
    slv_sw = TwoStreamSWRTE(grid_params; swbcs...)
    #------calling solvers
    metric_scaling = DA(one.(slv_sw.flux.flux_up))
    as = AtmosphericState(lon, lat, layerdata, p_lev, t_lev, t_sfc, vmr, cloud_state, nothing)

    n_bnd = 16
    flux_up_arr = zeros(FT,nlev,n_bnd)
    flux_dn_arr = zeros(FT,nlev,n_bnd)
    #βe = zeros(FT,nlev,256)
    solve_lw!(slv_lw, flux_up_arr, flux_dn_arr, as, lookup_lw, lookup_lw_cld, nothing, metric_scaling)
    solve_sw!(slv_sw, as, lookup_sw, lookup_sw_cld, nothing, metric_scaling)

    flux_net = slv_lw.flux.flux_net + slv_sw.flux.flux_net
    cp_d_ = FT(RRTMGP.Parameters.cp_d(param_set))
    grav_ = FT(RRTMGP.Parameters.grav(param_set))
    hr_lay = DA{FT}(undef, nlay, ncol)

    compute_gray_heating_rate!(device,hr_lay,p_lev,ncol,nlay,flux_net,cp_d_,grav_)
   
    if i==30
        println("here")
    end

    grid.states.T .+= -hr_lay[:,1] * dt 
    grid.states.θ .= theta_from_T(grid.states.T, grid.states.P, constants)
    grid.states.ρ .= ρ_ideal_gas(grid.states.P, grid.states.T, grid.states.qv, constants)
    
    ###Add in other updates
    
    flux_net_droplet = zeros(FT,nlay,n_bnd)
    for ibnd in 1:n_bnd
        
        totplnk = view(lookup_lw.planck.tot_planck, :, ibnd)
        (; t_planck) = lookup_lw.planck

        #extliq, ssaliq, _ = LookUpTables.getview_liqdata(lookup_lw_cld, ibnd)
      
        #re_vals = range(radliq_lwr,radliq_upr,lookup_lw_cld.dims[3])
        #ext_d = Optics.interp1d_equispaced(cld_r_eff_liq[glay,1],re_vals,extliq)
        #ssa_d = Optics.interp1d_equispaced(cld_r_eff_liq[glay,1],re_vals,ssaliq)

        #abs_d = (1.0-ssa_d) * ext_d
        
        @inbounds begin
            for glay in 1:nlay
                if cld_mask_lw[glay,1]
                    bb_flux = pi * Optics.interp1d_equispaced(t_lay[glay,1], t_planck, totplnk) 
                   
                    flux_net_droplet[glay,ibnd] = bb_flux - 0.5 * (flux_up_arr[glay,ibnd] + flux_dn_arr[glay+1,ibnd])
                    
                end
            end
        end
    end


    return flux_net_droplet
end