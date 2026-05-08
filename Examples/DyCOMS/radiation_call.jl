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


struct radiation_clima_data_helper{FT<:AbstractFloat} 
    context::ClimaComms.AbstractCommsContext
    device::ClimaComms.AbstractDevice
    DA::Type{<:AbstractArray}
    nlay::Int
    nlev::Int
    ncol::Int
    ngas::Int
    nband_lw::Int
    nband_sw::Int
    lookup_lw::LookUpLW
    lookup_lw_cld::LookUpCld
    lookup_sw::LookUpSW
    lookup_sw_cld::LookUpCld
    idx_gases::Dict{String,Int}
    idx_h2o::Int
    paramset::RRTMGPParameters
    grid_params::RRTMGPGridParams
    lat::Union{FT,Nothing}
    lon::Union{FT,Nothing}
    slv_lw::TwoStreamLWRTE
    slv_sw::TwoStreamSWRTE

    function radiation_clima_data_helper(nz::Int, FT::Type{<:AbstractFloat})
        nlay = nz + 10
        nlev = nlay + 1
        ncol = 1
        lookup_lw,lookup_lw_cld,lookup_sw,lookup_sw_cld,idx_gases = read_radiation_tables()
        context = ClimaComms.context()
        device = ClimaComms.device(context)
        DA = ClimaComms.array_type(device)
        param_set = RRTMGPParameters(9.80665, 0.028964, 0.018016, 8.3144598, 0.28571428571, 5.67e-8, 6.02214076e+23)
        grid_params = RRTMGPGridParams(FT; context, nlay, ncol)
        ngas = length(idx_gases)
        nband_lw = LookUpTables.get_n_bnd(lookup_lw)
        nband_sw = LookUpTables.get_n_bnd(lookup_sw)

        lat = nothing
        lon = nothing

    
        swbcs = (cos_zenith=[FT(-1.0)], toa_flux=[FT(lookup_sw.solar_src_tot)], sfc_alb_direct=[FT(0.06)], inc_flux_diffuse=nothing, sfc_alb_diffuse=[FT(0.06)])
        slv_sw = TwoStreamSWRTE(grid_params; swbcs...)
        
        # sfc_emis = fill(FT(0.98), 1)
        sfc_emis = fill(FT(0.98), 1)
        inc_flux = nothing
        slv_lw = TwoStreamLWRTE(grid_params; params = param_set, sfc_emis, inc_flux)
        
        return new{FT}(context,device,DA,nlay, nlev, ncol, ngas, nband_lw, nband_sw, lookup_lw, lookup_lw_cld, lookup_sw, 
        lookup_sw_cld, idx_gases, idx_gases["h2o"],
        param_set, grid_params, lat, lon, slv_lw, slv_sw)
    end
end


struct radiation_data{FT<:AbstractFloat}
    flux_up_lw::Array{FT,2}
    flux_dn_lw::Array{FT,2}
    flux_up_sw::Array{FT,2}
    flux_dn_sw::Array{FT,2}
    flux_net::Array{FT,2}
    layerdata::Array{FT, 2}
    metric_scaling::Array{FT, 2}
    flux_up_arr::Array{FT, 2}
    flux_dn_arr::Array{FT, 2}
    hr_lay::Array{FT, 2}
    atmospheric_state::AtmosphericState
    flux_net_droplet::Matrix{FT}
    cloud_heating_delta::Vector{FT}
    cond_rad_term::Vector{FT}
    CArad::radiation_clima_data_helper{FT}


    function radiation_data(::Type{FT},nsd,grid,constants) where FT<:AbstractFloat
        CArad = radiation_clima_data_helper(grid.nz, FT)
        nlay= CArad.nlay
        nlev = CArad.nlev
        ncol = CArad.ncol
        n_band_sw = CArad.nband_sw
        n_bnd = CArad.nband_lw
        levdim = (FT, nlev, ncol)
        laydim = (FT, nlay, ncol)

        #Profiles
        layerdata = similar(zeros(laydim...), 4, nlay)
        Nz = grid.nz
        θ, qv, P, P_faces = grid.states.θ, grid.states.qv, grid.states.P,grid.states.P_faces
        grid.states.T_tmp .= T_from_theta.(grid.states.θ,grid.states.P,constants)
        T = grid.states.T_tmp
        Tsfc = 292.0
        T1 = 0.5(T[1]+Tsfc)#2*T[1] - T[2] 
        P_extra = range(P_faces[end], 30000, 21)
        T_extra = T_from_theta.(θ[end],P_extra, constants) 
        layerdata[2,:] = [P; P_extra[2:2:end-1]]
        layerdata[3,:] = [T;T[end]; T_extra[4:2:end-1]]
        p_lev = reshape([P_faces; P_extra[3:2:end]], nlev, 1)
        t_lev = reshape([T1; (T[1:end-1] .+ T[2:end])./2; T_extra[1:2:end]], nlev, 1)

        #gas
        ngas = CArad.ngas
        idx_gases = CArad.idx_gases
        vmrat = zeros(FT, ngas, nlay, ncol)
        vmrat[idx_gases["o3"], :, 1] .= FT(0) #no ozone    
        vmrat[idx_gases["co2"], :, 1] .= FT(348e-6)
        vmrat[idx_gases["ch4"], :, 1] .= FT(1650e-9)
        vmrat[idx_gases["n2o"], :, 1] .= FT(306e-9)
        vmrat[idx_gases["n2"], :, 1] .= FT(0.7808)
        vmrat[idx_gases["o2"], :, 1] .= FT(0.2095)
        vmrat[idx_gases["co"], :, 1] .= FT(0)
        vmrat[idx_gases["h2o"], :, 1] .= [grid.states.qv;range(qv[end], 1e-8, 11)[2:end]] .* constants.ϵ
        vmr = Vmr(CArad.DA(vmrat))
        ice_rgh = 2
        #Clima Structs
        metric_scaling = CArad.DA(one.(CArad.slv_sw.flux.flux_up))
        cloud_state = CloudState(zeros(laydim...),zeros(laydim...),zeros(laydim...),zeros(laydim...),
                zeros(laydim...),falses(nlay,ncol),falses(nlay,ncol),MaxRandomOverlap(),ice_rgh)

        as =AtmosphericState(CArad.lon, CArad.lat, layerdata, p_lev, t_lev, Tsfc, vmr, cloud_state, nothing)        
        
        #output
        cloud_heating_delta = zeros(FT, Nz)
        cond_rad_term = zeros(FT, nsd)
        

        return new{FT}(zeros(levdim...),zeros(levdim...),zeros(levdim...),zeros(levdim...),zeros(levdim...),
        layerdata, metric_scaling, 
        zeros(FT,nlev,n_bnd),zeros(FT,nlev,n_bnd), zeros(laydim...), as, 
        zeros(FT, Nz, n_bnd), cloud_heating_delta, cond_rad_term, CArad)

        
    end 
end

Base.broadcastable(x::radiation_data) = Ref(x)
Base.broadcastable(x::radiation_clima_data_helper)= Ref(x)

function fill_full_atmosphere(grid::scm_eulerian_arrays,as::AtmosphericState,constants::Constants,CArad::radiation_clima_data_helper)::Nothing
    Nz = grid.nz
    grid.states.T_tmp .= T_from_theta.(grid.states.θ,grid.states.P,constants)
    T = grid.states.T_tmp
    T1 = 0.5(T[1]+292.0)# 2*T[1] - T[2]
    t_lay = @view as.layerdata[3, :, 1]
    t_lay[1:Nz] .= T
    as.t_lev[1, 1] = T1
    as.t_lev[2:Nz, 1] .= (@view(T[1:end-1]) .+ @view(T[2:end])) ./ 2
    # as.t_sfc         .= T1
    as.vmr.vmr[CArad.idx_h2o, 1:Nz, 1] .= grid.states.qv .* constants.ϵ

return nothing
end

function fill_liquid_water_diagnostics(grid, diagnosticsettings,CArad,as) 
    Nz = grid.nz
    radliq_lwr = CArad.lookup_lw_cld.bounds[1]
    radliq_upr = CArad.lookup_lw_cld.bounds[2]
    as.cloud_state.mask_lw .= false
    as.cloud_state.cld_frac .= 0.0

    as.cloud_state.cld_path_liq[1:Nz,1] .= grid.diagnostics.cloud_LWC * grid.dz * 1000
    as.cloud_state.cld_r_eff_liq[1:Nz,1] .= grid.diagnostics.cloud_effective_radius * 1.0e6
    as.cloud_state.cld_r_eff_liq[as.cloud_state.cld_r_eff_liq .< radliq_lwr] .= radliq_lwr
    as.cloud_state.cld_r_eff_liq[as.cloud_state.cld_r_eff_liq .> radliq_upr] .= radliq_upr
    as.cloud_state.cld_frac[as.cloud_state.cld_path_liq.>0.0] .= 1.0
    as.cloud_state.mask_lw[as.cloud_state.cld_frac.==1] .= true
    as.cloud_state.mask_sw .= as.cloud_state.mask_lw

    return nothing
end

function radiation_function!(grid,spatialsettings, diagnosticsettings, constants,raddata, dt, i)::Nothing
    
    #Some things to consider:
    #1. Do we want incoming downwelling LW?
    #2. Or, do we want to specify the atmosphere above the model top just for radiation?


    FT = eltype(raddata.flux_up_lw)
    CArad = raddata.CArad
    context = CArad.context
    device = CArad.device
    DA = CArad.DA
    param_set = CArad.paramset
    grid_params = CArad.grid_params
    as = raddata.atmospheric_state


    #Profiles
    fill_full_atmosphere(grid,as,constants,CArad)
    fill_liquid_water_diagnostics(grid, diagnosticsettings,CArad,as)

    nlay = CArad.nlay
    ncol = CArad.ncol
    n_bnd = CArad.nband_lw
    col_dry = reshape(view(as.layerdata, 1, :), nlay, ncol)
    p_lay = reshape(view(as.layerdata, 2, :), nlay, ncol)
    t_lay = reshape(view(as.layerdata, 3, :), nlay, ncol)
    rel_hum = reshape(view(as.layerdata, 4, :), nlay, ncol)
    vmr_h2o = reshape(view(as.vmr.vmr, CArad.idx_gases["h2o"], :, 1), nlay, ncol)

    compute_col_gas!(device, as.p_lev, col_dry, param_set, vmr_h2o, CArad.lat) # the example skips lat based gravity calculation
    compute_relative_humidity!(device, rel_hum, p_lay, t_lay, param_set, vmr_h2o) # compute relative humidity


    #βe = zeros(FT,nlev,256)
    solve_lw!(CArad.slv_lw, raddata.flux_up_arr, raddata.flux_dn_arr, as, CArad.lookup_lw, CArad.lookup_lw_cld, nothing, raddata.metric_scaling)
    # solve_sw!(CArad.slv_sw, as, CArad.lookup_sw, CArad.lookup_sw_cld, nothing,raddata.metric_scaling)

    raddata.flux_net .= CArad.slv_lw.flux.flux_net #.+ CArad.slv_sw.flux.flux_net
    cp_d_ = FT(RRTMGP.Parameters.cp_d(param_set))
    grav_ = FT(RRTMGP.Parameters.grav(param_set))

    compute_gray_heating_rate!(device,raddata.hr_lay,as.p_lev,ncol,nlay,raddata.flux_net,cp_d_,grav_)

    grid.states.T_tmp .+= dt .* @view raddata.hr_lay[1:nz,1] 
    raddata.cloud_heating_delta .+= dt .* @view raddata.hr_lay[1:nz,1]

    grid.states.θ .= theta_from_T(grid.states.T_tmp, grid.states.P, constants)
    grid.states.ρ .= ρ_ideal_gas(grid.states.P, grid.states.T_tmp, grid.states.qv, constants)
    
    ###Add in other updates
    
    raddata.flux_net_droplet .= 0#zeros(FT,nlay,n_bnd)
    for ibnd in 1:n_bnd
        
        totplnk = view(CArad.lookup_lw.planck.tot_planck, :, ibnd)
        (; t_planck) = CArad.lookup_lw.planck

        
        @inbounds begin
            for glay in 1:grid.nz
                if as.cloud_state.mask_lw[glay,1]
                    bb_flux = pi * Optics.interp1d_equispaced(t_lay[glay,1], t_planck, totplnk) # t_lay from the timestep before update
                   
                    # raddata.flux_net_droplet[glay,ibnd] = bb_flux - 0.5 * (raddata.flux_up_arr[glay,ibnd] + raddata.flux_dn_arr[glay+1,ibnd])
                    raddata.flux_net_droplet[glay,ibnd] = 4* bb_flux - (raddata.flux_up_arr[glay,ibnd] + raddata.flux_dn_arr[glay+1,ibnd])
                    
                end
            end
        end
    end


    return nothing#flux_net_droplet
end