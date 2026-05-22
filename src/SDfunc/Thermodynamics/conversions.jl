
# export calc_θ_dry, calc_ρ_dry_from_ρ, calc_T, calc_p
export theta_from_T, T_virtual, T_from_theta, ρ_ideal_gas, ρ_calc_θ,ρ_calc_θ!
export compute_ql_at_cell!,θl
#_θ, θl
export mixing_ratio, specific_humidity
export theta_from_T!, ρ_ideal_gas!, T_from_theta!

# function calc_θ_dry(constants,θl, q_vap)
#     FT = eltype(q_vap)
#     r_vap::FT = q_vap / (1 - q_vap)

#     return θl * (1 + r_vap / constants.ϵ)^(constants.Rd / constants.Cp_air)
# end

# function calc_ρ_dry_from_ρ(ρ, q_vap)
#     FT = eltype(ρ)
#     r_vap::FT = q_vap / (1 - q_vap)
#     return ρ / (1 + r_vap)
# end

# function calc_T(constants,θ_dry, ρ_dry)
#     FT = eltype(θ_dry)

#     R_d::FT = constants.Rd
#     cp_d::FT = constants.Cp_air
#     p_0::FT = constants.P0

#     return θ_dry * (ρ_dry * θ_dry / p_0 * R_d)^(R_d / cp_d / (1 - R_d / cp_d))
# end

# function calc_p(constants,ρ_dry, T, q_vap)
#     FT = eltype(T)
#     r_vap::FT = q_vap / (1 - q_vap)

#     R_d::FT = constants.Rd
#     R_v::FT = constants.Rv

#     return ρ_dry * (1 + r_vap) * (R_v / (1 / r_vap + 1) + R_d / (1 + r_vap)) * T
# end


#exner function, all of the translations
theta_from_T(T,P,constants) = T .* (constants.P0 ./ P).^(constants.Rd / constants.Cp_air)
T_virtual(T,q_vap) = T .* (1 .+ 0.61 .* q_vap)

θl(P,T,ql,constants) = (constants.P0 ./ P).^(constants.Rd / constants.Cp_air) .* (T .- constants.L .*ql ./constants.Cp_air) #should the last cp be for liquid?
# θl_θ(θ,ql,constants) = θ .- constants.L .*ql ./constants.Cp_air
ρ_ideal_gas(P,T,q_vap,constants) = P ./ (constants.Rd .* T_virtual(T,q_vap))
@inline ρ_calc_θ(P, θ, q_vap, constants) = ρ_ideal_gas(P,T_from_theta(θ,P,constants),q_vap,constants)


theta_from_T!(θ, T, P, constants) = @. θ = T * (constants.P0 / P)^(constants.Rd / constants.Cp_air)
ρ_calc_θ!(ρ, P, θ, q_vap, constants) =
    @. ρ = ρ_ideal_gas(P, T_from_theta(θ,P,constants), q_vap, constants)

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

@inline T_from_theta(θ,P,constants) = θ * (P / constants.P0)^(constants.Rd / constants.Cp_air)
T_from_theta!(T, θ, P, constants) = @. T = θ * (P / constants.P0)^(constants.Rd / constants.Cp_air)