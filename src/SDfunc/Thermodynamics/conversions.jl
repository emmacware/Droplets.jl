
export calc_θ_dry, calc_ρ_dry_from_ρ, calc_P_dry_from_P, calc_ρ_from_ρ_dry, calc_P_from_P_dry, calc_T#, calc_p
export theta_from_T, T_virtual, θ_virtual, θ_from_θv, T_from_theta, ρ_ideal_gas, ρ_calc_θ,ρ_calc_θ!
export compute_ql_at_cell!,compute_ql_all_cells!,θl,θ_from_θl
#_θ, θl
export mixing_ratio, specific_humidity
export theta_from_T!, ρ_ideal_gas!, T_from_theta!
export virtual_temp_coeff
  function calc_θ_dry(constants,θ, q_vap)
      r_vap = mixing_ratio(q_vap)
      return θ * (1 + r_vap / constants.ϵ)^(constants.Rd / constants.Cp_air)
  end
  function calc_θ_moist(constants,θd, q_vap)
      r_vap = mixing_ratio(q_vap)
      return θd * (1 + r_vap / constants.ϵ)^(-constants.Rd / constants.Cp_air)
  end

function calc_ρ_dry_from_ρ(ρ, q_vap)
    FT = eltype(ρ)
    r_vap::FT = q_vap / (1 - q_vap)
    return ρ / (1 + r_vap)
end


function calc_P_dry_from_P(P, q_vap, constants)
    FT = eltype(P)
    r_vap::FT = q_vap / (1 - q_vap)
    return P / (1 + r_vap / constants.ϵ)
end

function calc_ρ_from_ρ_dry(ρ_dry, q_vap)
    FT = eltype(ρ_dry)
    r_vap::FT = q_vap / (1 - q_vap)
    return ρ_dry * (1 + r_vap)
end

function calc_P_from_P_dry(P_dry, q_vap, constants)
    FT = eltype(P_dry)
    r_vap::FT = q_vap / (1 - q_vap)
    return P_dry * (1 + r_vap / constants.ϵ)
end

function calc_T(constants,θ, ρ_dry, q_vap)
    FT = eltype(θ)
    θ_dry::FT = calc_θ_dry(constants,θ, q_vap)

    R_d::FT = constants.Rd
    cp_d::FT = constants.Cp_air
    p_0::FT = constants.P0

    return θ_dry * (ρ_dry * θ_dry / p_0 * R_d)^(R_d / cp_d / (1 - R_d / cp_d))
end


# virtual T
# Virtual-temperature buoyancy coefficient Rv/Rd - 1 (0.61)
@inline virtual_temp_coeff(constants) = constants.Rv / constants.Rd - 1
T_virtual(T,q_vap,constants) = T .* (1 .+ virtual_temp_coeff(constants) .* q_vap)
θ_virtual(θ,q_vap,constants) = θ .* (1 .+ virtual_temp_coeff(constants) .* q_vap)
θ_from_θv(θv,q_vap,constants) = θv ./ (1 .+ virtual_temp_coeff(constants) .* q_vap)

@inline ρ_calc_θ(P, θ, q_vap, constants) = ρ_ideal_gas(P,T_from_theta(θ,P,q_vap,constants),q_vap,constants)
ρ_calc_θ!(ρ, P, θ, q_vap, constants) =
    @. ρ = ρ_ideal_gas(P, T_from_theta(θ,P,q_vap,constants), q_vap, constants)




function compute_ql_at_cell!(state, k, constants)
    r = state.droplets.grid_range[k]
    if isempty(r)
        state.ql_tmp[k] = zero(eltype(state.ρ))
        return
    end
    ids = @view state.droplets.I[r]
    s = zero(eltype(state.ρ))
    for j in ids
        s += state.droplets.X[j] * state.droplets.ξ[j]
    end
    state.ql_tmp[k] = s * constants.ρl / (state.ρ[k] * state.spatial.area_per_grid * state.spatial.z_grid_height)
end

# Single O(N) pass over droplets.cell_id instead of grid range (for multiple reasons)
function compute_ql_all_cells!(state, constants)
    ql = state.ql_tmp
    fill!(ql, zero(eltype(ql)))
    droplets = state.droplets
    nz = length(ql)
    @inbounds for i in droplets.I
        k = droplets.cell_id[i]
        (k < 1 || k > nz) && continue
        ql[k] += droplets.X[i] * droplets.ξ[i]
    end
    scale = constants.ρl / (state.spatial.area_per_grid * state.spatial.z_grid_height)
    @. ql = ql * scale / state.ρ
end


@inline mixing_ratio(q_vap) = q_vap / (1 - q_vap)
@inline specific_humidity(r) = r / (1 + r)

# -------------------------------------------------------
# θ/T Poisson-relation family. Dispatches on constants.dry_theta_convention:
#   false (standard): total pressure P, matches case-spec θ values directly.
#   true  (PySDM/KiD dry): dry partial pressure P_dry -- only self-consistent
#     when the state was initialized via the ρ_dry/P_dry (hydrostatic_pysdm)
#     path, since that's what calibrates what a given θ number means.
@inline R_m(q_vap, constants) = constants.Rd
@inline Cp_m(q_vap, constants) = constants.Cp_air

θl(P,T,ql,q_vap,constants) = constants.dry_theta_convention ?
    (constants.P0 ./ calc_P_dry_from_P.(P,q_vap,Ref(constants))).^(constants.Rd / constants.Cp_air) .* (T .- constants.L .*ql ./constants.Cp_air) :
    (constants.P0 ./ P).^(R_m(q_vap, constants) / Cp_m(q_vap, constants)) .* (T .- constants.L .*ql ./Cp_m.(q_vap, constants))

θ_from_θl(P,θl_val,ql,q_vap,constants) = constants.dry_theta_convention ?
    θl_val .+ (constants.P0 ./ calc_P_dry_from_P.(P,q_vap,Ref(constants))).^(constants.Rd / constants.Cp_air) .* constants.L .* ql ./ constants.Cp_air :
    θl_val .+ (constants.P0 ./ P).^(R_m(q_vap, constants) / Cp_m(q_vap, constants)) .* constants.L .* ql ./ Cp_m.(q_vap, constants)

theta_from_T(T,P,q_vap,constants) = constants.dry_theta_convention ?
    T .* (constants.P0 ./ calc_P_dry_from_P.(P,q_vap,Ref(constants))).^(constants.Rd / constants.Cp_air) :
    T .* (constants.P0 ./ P).^(R_m(q_vap, constants) / Cp_m(q_vap, constants))

theta_from_T!(θ, T, P, q_vap,constants) = constants.dry_theta_convention ?
    (@. θ = (T * (constants.P0 / calc_P_dry_from_P(P,q_vap,constants))^(constants.Rd / constants.Cp_air))) :
    (@. θ = (T * (constants.P0 / P)^(R_m(q_vap, constants) / Cp_m(q_vap, constants))))

@inline T_from_theta(θ,P,q_vap,constants) = constants.dry_theta_convention ?
    θ * (calc_P_dry_from_P(P,q_vap,constants) / constants.P0)^(constants.Rd / constants.Cp_air) :
    θ * (P / constants.P0)^(R_m(q_vap, constants) / Cp_m(q_vap, constants))

T_from_theta!(T, θ, P, q_vap,constants) = constants.dry_theta_convention ?
    (@. T = θ * (calc_P_dry_from_P(P,q_vap,constants) / constants.P0)^(constants.Rd / constants.Cp_air)) :
    (@. T = θ * (P / constants.P0)^(R_m(q_vap, constants) / Cp_m(q_vap, constants)))

@inline R_m_moist(q_vap, constants) = (constants.Rd + mixing_ratio(q_vap) * constants.Rv) / (1 + mixing_ratio(q_vap))

ρ_ideal_gas(P,T,q_vap,constants) = constants.dry_theta_convention ?
    P ./ (R_m_moist.(q_vap, constants) .* T) :
    P ./ (constants.Rd .* T_virtual(T,q_vap,constants))

ρ_ideal_gas!(ρ, P, T, q_vap, constants) = constants.dry_theta_convention ?
    (@. ρ = P / (R_m_moist(q_vap, constants) * T)) :
    (@. ρ = P / (constants.Rd * T_virtual(T,q_vap,constants)))

