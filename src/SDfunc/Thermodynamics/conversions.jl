
export calc_θ_dry, calc_ρ_dry_from_ρ, calc_T#, calc_p
export theta_from_T, T_virtual, T_from_theta, ρ_ideal_gas, ρ_calc_θ,ρ_calc_θ!
export compute_ql_at_cell!,θl,θ_from_θl
#_θ, θl
export mixing_ratio, specific_humidity
export theta_from_T!, ρ_ideal_gas!, T_from_theta!

function calc_θ_dry(constants,θ, q_vap)
    FT = eltype(q_vap)
    # r_vap::FT = q_vap / (1 - q_vap)

    return θ * (1 + q_vap / constants.ϵ)^(constants.Rd / constants.Cp_air)
end

function calc_θ_moist(constants,θd, q_vap)
    FT = eltype(q_vap)
    # r_vap::FT = q_vap / (1 - q_vap)
    return θd * (1 + q_vap / constants.ϵ)^(-constants.Rd / constants.Cp_air)
end

function calc_ρ_dry_from_ρ(ρ, q_vap)
    FT = eltype(ρ)
    r_vap::FT = q_vap / (1 - q_vap)
    return ρ / (1 + r_vap)
end

function calc_T(constants,θ, ρ_dry, q_vap)
    FT = eltype(θ)
    θ_dry::FT = calc_θ_dry(constants,θ, q_vap)

    R_d::FT = constants.Rd
    cp_d::FT = constants.Cp_air
    p_0::FT = constants.P0

    return θ_dry * (ρ_dry * θ_dry / p_0 * R_d)^(R_d / cp_d / (1 - R_d / cp_d))
end

# function calc_p(constants,ρ_dry, T, q_vap)
#     FT = eltype(T)
#     r_vap::FT = q_vap / (1 - q_vap)

#     R_d::FT = constants.Rd
#     R_v::FT = constants.Rv

#     return ρ_dry * (1 + r_vap) * (R_v / (1 / r_vap + 1) + R_d / (1 + r_vap)) * T
# end


#exner function, all of the translations
theta_from_T(T,P,q_vap,constants) = T .* (constants.P0 ./ P).^(R_m(q_vap, constants) / Cp_m(q_vap, constants))
T_virtual(T,q_vap) = T .* (1 .+ 0.61 .* q_vap)

θl(P,T,ql,q_vap,constants) = (constants.P0 ./ P).^(R_m(q_vap, constants) / Cp_m(q_vap, constants)) .* (T .- constants.L .*ql ./Cp_m.(q_vap, constants))
# Inverse of θl: recover θ given θl and the (possibly newly-redistributed) liquid content.
# θl = θ - Π⁻¹·L·ql/Cp_m, so θ = θl + Π⁻¹·L·ql/Cp_m.
θ_from_θl(P,θl_val,ql,q_vap,constants) = θl_val .+ (constants.P0 ./ P).^(R_m(q_vap, constants) / Cp_m(q_vap, constants)) .* constants.L .* ql ./ Cp_m.(q_vap, constants)
ρ_ideal_gas(P,T,q_vap,constants) = P ./ (constants.Rd .* T_virtual(T,q_vap))
@inline ρ_calc_θ(P, θ, q_vap, constants) = ρ_ideal_gas(P,T_from_theta(θ,P,q_vap,constants),q_vap,constants)


theta_from_T!(θ, T, P, q_vap,constants) = @. θ = (T * (constants.P0 / P)^(R_m(q_vap, constants) / Cp_m(q_vap, constants)))
ρ_calc_θ!(ρ, P, θ, q_vap, constants) =
    @. ρ = ρ_ideal_gas(P, T_from_theta(θ,P,q_vap,constants), q_vap, constants)

ρ_ideal_gas!(ρ, P, T, q_vap, constants) = @. ρ = P / (constants.Rd * T_virtual(T,q_vap))


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


@inline mixing_ratio(q_vap) = q_vap / (1 - q_vap)
@inline specific_humidity(r) = r / (1 + r)

@inline R_m(q_vap, constants) = constants.Rd #* (1 + mixing_ratio(q_vap) / constants.ϵ) / (1 + mixing_ratio(q_vap))
@inline Cp_m(q_vap, constants) = constants.Cp_air #+ mixing_ratio(q_vap) * constants.Cp_vapor) / (1 + mixing_ratio(q_vap))

@inline T_from_theta(θ,P,q_vap,constants) = θ * (P / constants.P0)^(R_m(q_vap, constants) / Cp_m(q_vap, constants))
T_from_theta!(T, θ, P, q_vap,constants) = @. T = θ * (P / constants.P0)^(R_m(q_vap, constants) / Cp_m(q_vap, constants))
