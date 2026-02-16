
export calc_θ_dry, calc_ρ_dry_from_ρ, calc_T, calc_p, theta_from_T, T_virtual, T_from_theta, ρ_ideal_gas



function calc_θ_dry(constants,θl, q_vap)
    FT = eltype(q_vap)
    r_vap::FT = q_vap / (1 - q_vap)

    return θl * (1 + r_vap * constants.ϵ)^(constants.Rd / constants.Cp_air)
end

function calc_ρ_dry_from_ρ(ρ, q_vap)
    FT = eltype(ρ)
    r_vap::FT = q_vap / (1 - q_vap)
    return ρ / (1 + r_vap)
end

function calc_T(constants,θ_dry, ρ_dry)
    FT = eltype(θ_dry)

    R_d::FT = constants.Rd
    cp_d::FT = constants.Cp_air
    p_0::FT = constants.P0

    return θ_dry * (ρ_dry * θ_dry / p_0 * R_d)^(R_d / cp_d / (1 - R_d / cp_d))
end

function calc_p(constants,ρ_dry, T, q_vap)
    FT = eltype(T)
    r_vap::FT = q_vap / (1 - q_vap)

    R_d::FT = constants.Rd
    R_v::FT = constants.Rv

    return ρ_dry * (1 + r_vap) * (R_v / (1 / r_vap + 1) + R_d / (1 + r_vap)) * T
end


#exner function, all of the translations
theta_from_T(T,P,constants) = T .* (constants.P0 ./ P).^(constants.Rd / constants.Cp_air)
T_virtual(T,q_vap) = T .* (1 .+ 0.61 .* q_vap)
T_from_theta(θ,P,constants) = θ .* (P ./ constants.P0).^(constants.Rd / constants.Cp_air)
θl(P,T,ql,constants) = (constants.P0 ./ P).^(constants.Rd / constants.Cp_air) .* (T .- constants.L .*ql ./constants.Cp_air) #should the last cp be for liquid?
ρ_ideal_gas(P,T,q_vap,constants) = P ./ (constants.Rd .* T_virtual(T,q_vap))