using RRTMGP
using NCDatasets
using ClimaComms
using RRTMGP: RRTMGPGridParams
using RRTMGP.VolumeMixingRatios
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.RTESolver
using RRTMGP.GrayAtmosphere
import RRTMGP.Parameters.RRTMGPParameters

#these inputs have the same names as the scm.jl file so if you run that first you can
#call this function without it needing to be part of the timestepper

function read_radiation_tables()

    FT = get(ARGS, 1, Float64) == "Float32" ? Float32 : Float64
    context = ClimaComms.context()
    device = ClimaComms.device(context)
    DA = ClimaComms.array_type(device)

    # lw_file = joinpath(@__DIR__, "../../data/rrtmgp-data-lw-g256-2018-12-04.nc")
    lw_file = joinpath(@__DIR__, "../../data/rrtmgp-gas-lw-g256.nc")
    sw_file = joinpath(@__DIR__, "../../data/rrtmgp-data-sw-g224-2018-12-04.nc")
    # lw_cld_file = joinpath(@__DIR__, "../../data/rrtmgp-clouds-lw-g256.nc")
    lw_cld_file = joinpath(@__DIR__, "../../data/rrtmgp-clouds-lw-bnd.nc")
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


struct radiation_clima_data_helper{FT<:AbstractFloat, C, D, TDA, LLW, LLWC, LSW, LSWC, PS, GP, SLW, SSW}
    context::C
    device::D
    DA::TDA
    nlay::Int
    nlev::Int
    ncol::Int
    ngas::Int
    nband_lw::Int
    nband_sw::Int
    lookup_lw::LLW
    lookup_lw_cld::LLWC
    lookup_sw::LSW
    lookup_sw_cld::LSWC
    idx_gases::Dict{String,Int}
    idx_h2o::Int
    paramset::PS
    grid_params::GP
    lat::Union{FT,Nothing}
    lon::Union{FT,Nothing}
    slv_lw::SLW
    slv_sw::SSW

    function radiation_clima_data_helper(nz::Int, FT::Type{<:AbstractFloat})
        nlay = nz + 10
        nlev = nlay + 1
        ncol = 1
        lookup_lw,lookup_lw_cld,lookup_sw,lookup_sw_cld,idx_gases = read_radiation_tables()
        # This driver always runs ncol = 1 (single-column model); multithreading a
        # loop of length 1 has nothing to parallelize and only adds per-g-point
        # `Threads.@threads` task-spawn overhead (and allocations) to every solve.
        context = ClimaComms.context(ClimaComms.CPUSingleThreaded())
        device = ClimaComms.device(context)
        DA = ClimaComms.array_type(device)
        param_set = RRTMGPParameters(9.80665, 0.028964, 0.018016, 8.3144598, 0.28571428571, 5.67e-8, 6.02214076e+23)
        grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol)
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
        slv_lw0 = TwoStreamLWRTE(grid_params; params = param_set, sfc_emis, inc_flux)
        # Attach a per-band flux buffer so `solve_lw!` retains band-resolved
        # up/down fluxes (used below for the cloud radiative heating term),
        # matching the shape `(nlev, ncol, n_bnd)`.
        band_flux_lw = FluxBand(grid_params, nband_lw)
        slv_lw = TwoStreamLWRTE(
            slv_lw0.context, slv_lw0.op, slv_lw0.src, slv_lw0.bcs,
            slv_lw0.fluxb, slv_lw0.flux, band_flux_lw, slv_lw0.state_cache,
        )

        return new{
            FT, typeof(context), typeof(device), typeof(DA),
            typeof(lookup_lw), typeof(lookup_lw_cld), typeof(lookup_sw), typeof(lookup_sw_cld),
            typeof(param_set), typeof(grid_params), typeof(slv_lw), typeof(slv_sw),
        }(context,device,DA,nlay, nlev, ncol, ngas, nband_lw, nband_sw, lookup_lw, lookup_lw_cld, lookup_sw,
        lookup_sw_cld, idx_gases, idx_gases["h2o"],
        param_set, grid_params, lat, lon, slv_lw, slv_sw)
    end
end


struct radiation_data{FT<:AbstractFloat, AS, CR}
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
    atmospheric_state::AS
    flux_net_droplet::Matrix{FT}
    cloud_heating_delta::Vector{FT}
    cond_rad_term::Vector{FT}
    CArad::CR


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
        T_from_theta!(grid.states.T_tmp,grid.states.θ,grid.states.P,grid.states.qv,constants)
        T = grid.states.T_tmp
        Tsfc = 292.0
        T1 = 0.5*(T[1]+Tsfc)#2*T[1] - T[2] 
        P_extra = range(P_faces[end], 30000, 21)
        T_extra = T_from_theta.(θ[end],P_extra,qv[end], constants) 
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
        vmrat[idx_gases["h2o"], :, 1] .= [mixing_ratio.(qv);range(mixing_ratio.(qv)[end], mixing_ratio(1e-8), 11)[2:end]] ./ constants.ϵ
        vmr = Vmr(CArad.DA(vmrat))
        ice_rgh = 2
        #Clima Structs
        # `apply_metric_scaling!` expects a vertical-first `(nlev, ncol)` array,
        # regardless of the internal (column-first) layout of the flux compute buffers.
        metric_scaling = CArad.DA(ones(FT, nlev, ncol))
        cloud_state = CloudState(zeros(laydim...),zeros(laydim...),zeros(laydim...),zeros(laydim...),
                zeros(laydim...),falses(nlay,ncol),falses(nlay,ncol),MaxRandomOverlap(),ice_rgh)

        as =AtmosphericState(CArad.lon, CArad.lat, layerdata, p_lev, t_lev, Tsfc, vmr, cloud_state, nothing)        
        
        #output
        cloud_heating_delta = zeros(FT, Nz)
        cond_rad_term = zeros(FT, nsd)
        

        return new{FT, typeof(as), typeof(CArad)}(zeros(levdim...),zeros(levdim...),zeros(levdim...),zeros(levdim...),zeros(levdim...),
        layerdata, metric_scaling,
        zeros(FT,nlev,n_bnd),zeros(FT,nlev,n_bnd), zeros(laydim...), as,
        zeros(FT, Nz, n_bnd), cloud_heating_delta, cond_rad_term, CArad)

        
    end 
end

Base.broadcastable(x::radiation_data) = Ref(x)
Base.broadcastable(x::radiation_clima_data_helper)= Ref(x)

function fill_full_atmosphere(grid::scm_eulerian_arrays,as::AtmosphericState,constants::Constants,CArad::radiation_clima_data_helper)::Nothing
    Nz = grid.nz
    T_from_theta!(grid.states.T_tmp,grid.states.θ,grid.states.P,grid.states.qv,constants)
    T = grid.states.T_tmp
    T1 = 0.5*(T[1]+292.0)# 2*T[1] - T[2]
    # T1 = 2*T[1] - T[2]
    t_lay = @view as.layerdata[3, :, 1]
    t_lay[1:Nz] .= T
    as.t_lev[1, 1] = T1
    # as.t_lev[2:Nz, 1] .= (@view(T[1:end-1]) .+ @view(T[2:end])) ./ 2
    @views as.t_lev[2:Nz,1] .= (T[1:end-1] .+ T[2:end]) .* 0.5
    # as.t_sfc         .= T1
    as.vmr.vmr[CArad.idx_h2o, 1:Nz, 1] .=  mixing_ratio.(grid.states.qv) ./ constants.ϵ

    return nothing
end

function fill_liquid_water_diagnostics(grid::scm_eulerian_arrays{FT}, diagnosticsettings::diagnostic_settings{FT},CArad::radiation_clima_data_helper,as::AtmosphericState)::Nothing where FT<:AbstractFloat
    Nz = grid.nz

    as.cloud_state.cld_path_liq[1:Nz,1] .= grid.diagnostics.cloud_LWC .* grid.dz .* 1000
    as.cloud_state.cld_r_eff_liq[1:Nz,1] .= grid.diagnostics.cloud_effective_radius .* 1.0e6
    clamp!(as.cloud_state.cld_r_eff_liq, CArad.lookup_lw_cld.bounds[1], CArad.lookup_lw_cld.bounds[2])
    
    @. as.cloud_state.cld_frac = ifelse(as.cloud_state.cld_path_liq > 0.0, one(FT), zero(FT))
    @. as.cloud_state.mask_lw = as.cloud_state.cld_frac == 1
    as.cloud_state.mask_sw .= as.cloud_state.mask_lw

    return nothing
end

# Function barrier for the ibnd/glay loop in compute_radiative_fluxes!: `planck` and
# `mask_lw` arrive from callers holding them as abstractly-typed fields (see the note
# at the call site below), so this function's own parameters get freshly specialized
# on their actual concrete runtime types instead of inheriting that instability.
function compute_flux_net_droplet!(raddata, planck, mask_lw, t_lay, n_bnd::Int, nz::Int)
    t_planck = planck.t_planck
    tot_planck = planck.tot_planck
    raddata.flux_net_droplet .= 0
    @inbounds for ibnd in 1:n_bnd
        totplnk = view(tot_planck, :, ibnd)
        for glay in 1:nz
            mask_lw[glay, 1] || continue
            bb_flux = π * Optics.interp1d_equispaced(t_lay[glay, 1], t_planck, totplnk)
            raddata.flux_net_droplet[glay, ibnd] = bb_flux - 0.25 * (raddata.flux_up_arr[glay-1, ibnd] + raddata.flux_up_arr[glay, ibnd] + raddata.flux_dn_arr[glay, ibnd] + raddata.flux_dn_arr[glay+1, ibnd])
        end
    end
    return nothing
end

function compute_radiative_fluxes!(grid::scm_eulerian_arrays{FT},spatialsettings::spatial_settings_1d{FT}, diagnosticsettings::diagnostic_settings{FT}, constants::Constants{FT},raddata::radiation_data)::Nothing where FT<:AbstractFloat

    #Some things to consider:
    #1. Do we want incoming downwelling LW?
    #2. Or, do we want to specify the atmosphere above the model top just for radiation?
    raddata.flux_up_arr .= 0
    raddata.flux_dn_arr .= 0
    raddata.flux_net .= 0
    # raddata.hr_lay .= 0
    # raddata.flux_net_droplet .= 0
    # raddata.cond_rad_term .= 0
    raddata.flux_up_lw .= 0
    raddata.flux_dn_lw .= 0
    raddata.flux_up_sw .= 0
    raddata.flux_dn_sw .= 0



    # FT = eltype(raddata.flux_up_lw)
    CArad = raddata.CArad
    context = CArad.context
    device = CArad.device
    DA = CArad.DA
    param_set = CArad.paramset
    grid_params = CArad.grid_params
    as = raddata.atmospheric_state
    nz = grid.nz


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
    solve_lw!(CArad.slv_lw, as, CArad.lookup_lw, CArad.lookup_lw_cld, nothing, raddata.metric_scaling)
    # solve_sw!(CArad.slv_sw, as, CArad.lookup_sw, CArad.lookup_sw_cld, nothing,raddata.metric_scaling)

    # `solve_lw!` fills CArad.slv_lw.band_flux (nlev, ncol, n_bnd); ncol == 1 here.
    raddata.flux_up_arr .= view(CArad.slv_lw.band_flux.flux_up, :, 1, :)
    raddata.flux_dn_arr .= view(CArad.slv_lw.band_flux.flux_dn, :, 1, :)

    # `flux_net` is stored column-first (ncol, nlev); raddata expects vertical-first (nlev, ncol).
    # An explicit loop (rather than `transpose(...)` broadcast) avoids a few bytes of
    # per-call allocation from the lazy-transpose wrapper on Julia <= 1.11.
    flux_net_lw = CArad.slv_lw.flux.flux_net
    @inbounds for gcol in 1:ncol, ilev in 1:size(raddata.flux_net, 1)
        raddata.flux_net[ilev, gcol] = flux_net_lw[gcol, ilev]
    end
    #.+ CArad.slv_sw.flux.flux_net
    cp_d_ = constants.Cp_air#FT(RRTMGP.Parameters.cp_d(param_set))
    grav_ = constants.gconst#FT(RRTMGP.Parameters.grav(param_set))

    compute_gray_heating_rate!(device,raddata.hr_lay,as.p_lev,ncol,nlay,raddata.flux_net,cp_d_,grav_)

    ##Add in other updates
    # CArad.lookup_lw and raddata.atmospheric_state are both declared with their
    # abstract supertypes (no type params) in radiation_clima_data_helper/radiation_data,
    # so accessing them -- and anything chained off them, e.g. t_planck below -- is
    # type-unstable. Rather than reparametrizing those structs (which would ripple
    # through every signature that names them), push the hot ibnd/glay loop through a
    # function barrier: compute_flux_net_droplet! gets called once with these
    # abstractly-typed values, but Julia specializes its *body* on their actual
    # concrete runtime types, so the loop inside it is fully concrete.
    compute_flux_net_droplet!(raddata, CArad.lookup_lw.planck, as.cloud_state.mask_lw, t_lay, n_bnd, grid.nz)

    return nothing
end

# Applies the last-computed radiative heating rate (raddata.hr_lay) to T/θ/ρ for one
# timestep. Runs every step regardless of whether compute_radiative_fluxes! was called
function apply_radiative_heating!(grid::scm_eulerian_arrays{FT}, raddata::radiation_data, constants::Constants{FT}, dt::FT)::Nothing where FT<:AbstractFloat
    nz = grid.nz
    hr = @view raddata.hr_lay[1:nz,1]

    grid.states.T_tmp .+= dt .* hr
    raddata.cloud_heating_delta .+= dt .* hr

    theta_from_T!(grid.states.θ, grid.states.T_tmp, grid.states.P,grid.states.qv, constants)
    ρ_ideal_gas!(grid.states.ρ, grid.states.P, grid.states.T_tmp, grid.states.qv, constants)

    return nothing
end

# Extends Droplets.radiation_function so that Droplets isnt dependent on RRTMGP
function Droplets.radiation_function!(::DynON,grid::scm_eulerian_arrays{FT},spatialsettings::spatial_settings_1d{FT}, diagnosticsettings::diagnostic_settings{FT}, constants::Constants{FT},raddata::radiation_data, dt::FT, i::Int, n_rad::Int=1)::Nothing where FT<:AbstractFloat
    if (i) % n_rad == 0 || i == 1
        compute_radiative_fluxes!(grid, spatialsettings, diagnosticsettings, constants, raddata)
    end

    apply_radiative_heating!(grid, raddata, constants, dt)

    return nothing
end

