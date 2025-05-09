include(joinpath(pkgdir(CM), "parcel", "ParcelParameters.jl"))

# """
#     If negative, clip the tracers to zero when computing tendencies
# """
clip!(x, lim) = max(x, lim)
clip!(x) = clip!(x, eltype(x)(0))

"""
    ODE problem definitions, using superdroplets
"""
function parcel_model_sd(dY, Y_CompVec, p, t)
    # Numerical precision used in the simulation
    Y = Y_CompVec.IC
    FT = eltype(Y)
    # Simulation parameters
    (; prescribed_thermodynamics, t_profile, T_profile, P_profile) = p
    (; wps, tps, r_nuc, w) = p
    (; liq_distr, ice_distr) = p
    (; aero_act_params, dep_params, imm_params, hom_params) = p
    (; ce_params, ds_params) = p
    # Y values stored in a named tuple for ease of use
    state = (
        Sₗ = Y[1],
        p_air = Y[2],
        T = Y[3],
        qᵥ = clip!(Y[4]),
        qₗ = clip!(Y[5]),
        qᵢ = clip!(Y[6]),
        Nₐ = clip!(Y[7]),
        Nₗ = clip!(Y[8]),
        Nᵢ = clip!(Y[9]),
        ln_INPC = Y[10],   # needed only in stochastic Frostenberg
        t = t,
    )

    # Constants
    Rᵥ = TD.Parameters.R_v(tps)
    grav = TD.Parameters.grav(tps)
    ρᵢ = wps.ρi
    ρₗ = wps.ρw

    # Get the state values
    (; Sₗ, p_air, T, qᵥ, qₗ, qᵢ, Nₗ, Nᵢ, t) = state
    # Get thermodynamic parameters, phase partition and create thermo state.
    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    ts = TD.PhaseNonEquil_pTq(tps, p_air, T, q)

    # Constants and variables that depend on the moisture content
    R_air = TD.gas_constant_air(tps, q)
    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)
    L_fus = TD.latent_heat_fusion(tps, T)
    L_vap = TD.latent_heat_vapor(tps, T)
    ρ_air = TD.air_density(tps, ts)

    # Adiabatic parcel coefficients
    a1 = L_vap * grav / cp_air / T^2 / Rᵥ - grav / R_air / T
    a2 = 1 / qᵥ
    a3 = L_vap^2 / Rᵥ / T^2 / cp_air
    a4 = L_vap * L_subl / Rᵥ / T^2 / cp_air
    a5 = L_vap * L_fus / Rᵥ / cp_air / (T^2)

    # Mean radius, area and volume of liquid droplets and ice crystals
    PSD_liq = distribution_moments(liq_distr, qₗ, Nₗ, ρₗ, ρ_air)
    PSD_ice = distribution_moments(ice_distr, qᵢ, Nᵢ, ρᵢ, ρ_air)

    # Aerosol activation
    dNₗ_dt_act = aerosol_activation(aero_act_params, state)
    dqₗ_dt_act = dNₗ_dt_act * 4 * FT(π) / 3 * r_nuc^3 * ρₗ / ρ_air

    # Deposition ice nucleation
    # (All deposition parameterizations assume monodisperse aerosol size distr)
    dNᵢ_dt_dep = deposition_nucleation(dep_params, state, dY.IC)
    dqᵢ_dt_dep = dNᵢ_dt_dep * 4 / 3 * FT(π) * r_nuc^3 * ρᵢ / ρ_air

    # Heterogeneous ice nucleation
    #@info(imm_params)
    dln_INPC_imm = INPC_model(imm_params, state)
    dNᵢ_dt_imm = immersion_freezing(imm_params, PSD_liq, state)
    dqᵢ_dt_imm = dNᵢ_dt_imm * PSD_liq.V * ρᵢ / ρ_air

    # Homogeneous ice nucleation
    dNᵢ_dt_hom = homogeneous_freezing(hom_params, PSD_liq, state)
    dqᵢ_dt_hom = dNᵢ_dt_hom * PSD_liq.V * ρᵢ / ρ_air

    # Condensation/evaporation
    R = volume_to_radius.(Y_CompVec.X)
    dX = dXkohler_function_of_radius_activated.(R,T,Sₗ,t)

    if ce_params == "Superdroplets"
        dqₗ_dt_ce = sum(dX.*Y_CompVec.ξ.* ρₗ / ρ_air) #/per volume?
    else
        dqₗ_dt_ce = condensation(ce_params, PSD_liq, state, ρ_air)
    end

    # Deposition/sublimation
    dqᵢ_dt_ds = deposition(ds_params, PSD_ice, state, ρ_air)

    # number concentration and ...
    dNᵢ_dt = dNᵢ_dt_dep + dNᵢ_dt_imm + dNᵢ_dt_hom
    dNₐ_dt = -dNᵢ_dt_dep - dNₗ_dt_act
    dNₗ_dt = dNₗ_dt_act - dNᵢ_dt_imm - dNᵢ_dt_hom
    # ... water mass budget
    dqₗ_dt_v2l = dqₗ_dt_ce + dqₗ_dt_act
    dqᵢ_dt_l2i = dqᵢ_dt_imm + dqᵢ_dt_hom
    dqᵢ_dt_v2i = dqᵢ_dt_dep + dqᵢ_dt_ds

    # Update the tendecies
    dqᵢ_dt = dqᵢ_dt_v2i + dqᵢ_dt_l2i
    dqₗ_dt = dqₗ_dt_v2l - dqᵢ_dt_l2i
    dqᵥ_dt = -dqₗ_dt_v2l - dqᵢ_dt_v2i

    dSₗ_dt =
        a1 * w * Sₗ - (a2 + a3) * Sₗ * dqₗ_dt_v2l -
        (a2 + a4) * Sₗ * dqᵢ_dt_v2i - a5 * Sₗ * dqᵢ_dt_l2i

    dp_air_dt =
        prescribed_thermodynamics ? AIDA_rate(t, t_profile, P_profile) :
        -p_air * grav / R_air / T * w

    dT_dt =
        prescribed_thermodynamics ? AIDA_rate(t, t_profile, T_profile) :
        -grav / cp_air * w +
        L_vap / cp_air * dqₗ_dt_v2l +
        L_fus / cp_air * dqᵢ_dt_l2i +
        L_subl / cp_air * dqᵢ_dt_v2i

    # Set tendencies
    dY[1] = dSₗ_dt      # saturation ratio over liquid water
    dY[2] = dp_air_dt   # pressure
    dY[3] = dT_dt       # temperature
    dY[4] = dqᵥ_dt      # vapor specific humidity
    dY[5] = dqₗ_dt      # liquid water specific humidity
    dY[6] = dqᵢ_dt      # ice specific humidity
    dY[7] = dNₐ_dt      # number concentration of interstitial aerosol
    dY[8] = dNₗ_dt      # mumber concentration of droplets
    dY[9] = dNᵢ_dt      # number concentration of activated particles
    dY[10] = dln_INPC_imm
    dY.X = dX
    dY.ξ .= 0
end

"""
    run_parcel_sd(IC, t_0, t_end, pp)

Returns the solution of an ODE probelm defined by the parcel model, with Superdroplet tracers.
"""

function run_parcel_sd(d_Y, t_0, t_end, pp)

    FT = eltype(d_Y.IC)

    info = "\nSize distribution (liq): $(pp.liq_size_distribution)\n"
    if pp.liq_size_distribution == "Monodisperse"
        liq_distr = Monodisperse{FT}()
    elseif pp.liq_size_distribution == "Gamma"
        liq_distr = Gamma{FT}()
    else
        throw("Unrecognized size distribution")
    end

    info *= "Size distribution (ice): $(pp.ice_size_distribution)\n"
    if pp.ice_size_distribution == "Monodisperse"
        ice_distr = Monodisperse{FT}()
    elseif pp.ice_size_distribution == "Gamma"
        ice_distr = Gamma{FT}()
    else
        throw("Unrecognized size distribution")
    end

    info *= "Aerosol: $(chop( string(typeof(pp.aerosol)), head = 29, tail = 9))\n"

    info *= "Aerosol Activation: $(pp.aerosol_act)\n"
    if pp.aerosol_act == "None"
        aero_act_params = Empty{FT}()
    elseif pp.aerosol_act == "AeroAct"
        aero_act_params =
            AeroAct{FT}(pp.aap, pp.aerosol, pp.aero_σ_g, pp.r_nuc, pp.const_dt)
    else
        throw("Unrecognized aerosol activation mode")
    end

    info *= "Deposition: $(pp.deposition)\n"
    if pp.deposition == "None"
        dep_params = Empty{FT}()
    elseif pp.deposition == "MohlerAF"
        dep_params = MohlerAF{FT}(pp.ips, pp.aerosol, pp.tps, pp.const_dt)
    elseif pp.deposition == "MohlerRate"
        dep_params = MohlerRate{FT}(pp.ips, pp.aerosol, pp.tps, pp.const_dt)
    elseif pp.deposition == "ABDINM"
        dep_params = ABDINM{FT}(pp.tps, pp.aerosol, pp.r_nuc, pp.const_dt)
    elseif pp.deposition == "P3_dep"
        dep_params = P3_dep{FT}(pp.ips, pp.const_dt)
    else
        throw("Unrecognized deposition mode")
    end

    info *= "Heterogeneous: $(pp.heterogeneous)\n"
    if pp.heterogeneous == "None"
        imm_params = Empty{FT}()
    elseif pp.heterogeneous == "ABIFM"
        imm_params = ABIFM{FT}(pp.tps, pp.aerosol, pp.A_aer, pp.const_dt)
    elseif pp.heterogeneous == "P3_het"
        imm_params = P3_het{FT}(pp.ips, pp.const_dt)
    elseif pp.heterogeneous == "Frostenberg_random"
        imm_params =
            Frostenberg_random{FT}(pp.ip, pp.sampling_interval, pp.const_dt)
    elseif pp.heterogeneous == "Frostenberg_mean"
        imm_params = Frostenberg_mean{FT}(pp.ip, pp.const_dt)
    elseif pp.heterogeneous == "Frostenberg_stochastic"
        imm_params = Frostenberg_stochastic{FT}(pp.ip, pp.γ, pp.const_dt)
    else
        throw("Unrecognized heterogeneous mode")
    end

    info *= "Homogeneous: $(pp.homogeneous)\n"
    if pp.homogeneous == "None"
        hom_params = Empty{FT}()
    elseif pp.homogeneous == "ABHOM"
        hom_params = ABHOM{FT}(pp.tps, pp.ips, pp.const_dt)
    elseif pp.homogeneous == "P3_hom"
        hom_params = P3_hom{FT}(pp.const_dt)
    else
        throw("Unrecognized homogeneous mode")
    end

    info *= "Condensation growth: $(pp.condensation_growth)\n"
    if pp.condensation_growth == "None"
        ce_params = Empty{FT}()
    elseif pp.condensation_growth == "Condensation"
        ce_params = CondParams{FT}(pp.aps, pp.tps, pp.const_dt)
    elseif pp.condensation_growth == "Superdroplets"
        ce_params = "Superdroplets"
    else
        throw("Unrecognized condensation growth mode")
    end

    info *= "Deposition growth: $(pp.deposition_growth)\n"
    if pp.deposition_growth == "None"
        ds_params = Empty{FT}()
    elseif pp.deposition_growth == "Deposition"
        ds_params = DepParams{FT}(pp.aps, pp.tps)
    else
        throw("Unrecognized deposition growth mode")
    end
    # To see parcel parameters in output, uncomment
    # @info info

    # Parameters for the ODE solver
    p = (
        prescribed_thermodynamics = pp.prescribed_thermodynamics,
        t_profile = pp.t_profile,
        T_profile = pp.T_profile,
        P_profile = pp.P_profile,
        liq_distr = liq_distr,
        ice_distr = ice_distr,
        aero_act_params = aero_act_params,
        dep_params = dep_params,
        imm_params = imm_params,
        hom_params = hom_params,
        ce_params = ce_params,
        ds_params = ds_params,
        wps = pp.wps,
        tps = pp.tps,
        aps = pp.aps,
        aap = pp.aap,
        r_nuc = pp.r_nuc,
        w = pp.w,
    )

    problem = ODE.ODEProblem(parcel_model_sd, d_Y, (FT(t_0), FT(t_end)), p)
    return ODE.solve(
        problem,
        ODE.Euler(),
        dt = pp.const_dt,#/100,
        reltol = 100 * eps(FT),
        abstol = 100 * eps(FT),
    )
end
