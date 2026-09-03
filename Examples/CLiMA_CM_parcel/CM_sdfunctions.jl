include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

"""
    SDCondensation{FT} <: CMP.ParametersType

`ce_params` type dispatched on by CM's real `condensation` (from `ParcelTendencies.jl`).
Holds the superdroplet volumes so `condensation` can compute the per-droplet Kohler growth
rate and aggregate it into the bulk qₗ tendency, instead of the bulk `CondParams` scheme.

`X` is re-synced from the ODE state each RHS call; `ξ` is fixed for the run (no coalescence)
dX is used for the wrapper parcel_model_sd
"""
mutable struct SDCondensation{FT} <: CMP.ParametersType
    X::Vector{FT}
    ξ::Vector{FT}
    dX::Vector{FT}
    ρₗ::FT
    dt::FT
end

function condensation(params::SDCondensation, PSD_liq, state, ρ_air)
    R = volume_to_radius.(params.X)
    params.dX .= dXkohler_function_of_radius_activated.(R, state.T, state.Sₗ, constants, params.dt)
    return sum(params.dX .* params.ξ) * params.ρₗ / ρ_air
end

"""
    parcel_model_sd(dY, Y, p, t)

Wrapper around CM's parcel example, unmodified parcel_model. Because Y = ComponentVector(IC=IC, ξ=..., X=...)
indexes CM's IC as the first 10 elements of the flat backing array, the positional Y[1]`...`Y[10] 
of the originial function remains correct. This wrapper only has to sync the current superdroplet volumes into
ce_params before the call (so the condensation function works), and write the resulting per-droplet growth 
rate into `dY.X` after.

"""
function parcel_model_sd(dY, Y, p, t)
    p.ce_params.X .= Y.X
    parcel_model(dY, Y, p, t)
    dY.X .= p.ce_params.dX
    dY.ξ .= 0
    return dY
end

"""
run_parcel_sd(d_Y, t_0, t_end, pp::parcel_params{FT})

only condensation is coupled; all other parcel_params process switches must be "None"
"""
function run_parcel_sd(d_Y, t_0, t_end, pp::parcel_params{FT}) where {FT}
    @assert pp.aerosol_act == "None" && pp.deposition == "None" &&
            pp.heterogeneous == "None" && pp.homogeneous == "None" &&
            pp.deposition_growth == "None" "run_parcel_sd only couples the condensation term; all other parcel_params process switches must be \"None\""

    # Droplets uses Exponential not Gamma, but CM doesnt have Exponential(), so we use Gamma() here .. doesnt do anything since we use Droplets condensation function and no ice
    liq_distr = pp.liq_size_distribution == "Monodisperse" ? Monodisperse() :
                pp.liq_size_distribution == "Exponential" ? Gamma() : 
                throw("Unrecognized liquid size distribution")

    ce_params = SDCondensation{FT}(copy(d_Y.X), FT.(d_Y.ξ), similar(d_Y.X), pp.wps.ρw, FT(pp.const_dt))

    p = (
        prescribed_thermodynamics = pp.prescribed_thermodynamics,
        t_profile = pp.t_profile,
        T_profile = pp.T_profile,
        P_profile = pp.P_profile,
        liq_distr = liq_distr,
        ice_distr = Monodisperse(), # ice is not functional but parcel_model still evaluates PSD_ice unconditionally
        aero_act_params = Empty(),
        dep_params = Empty(),
        imm_params = Empty(),
        hom_params = Empty(),
        ce_params = ce_params,
        ds_params = Empty(),
        wps = pp.wps,
        tps = pp.tps,
        r_nuc = pp.r_nuc,
        w = pp.w,
    )

    problem = ODE.ODEProblem(parcel_model_sd, d_Y, (FT(t_0), FT(t_end)), p)
    return ODE.solve(
        problem,
        ODE.Euler(),
        dt = pp.const_dt,
        reltol = 100 * eps(FT),
        abstol = 100 * eps(FT),
    )
end
