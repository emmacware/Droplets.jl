
function drho_dz(ρ, ode_settings, z)

    FT = eltype(ρ)
    constants = ode_settings

    θ_z::FT = θl(z)
    q_vap::FT = qt(z)

    # constants
    g::FT = constants.gconst
    R_d::FT = constants.Rd
    R_v::FT = constants.Rv
    cp_d::FT = constants.Cp_air
    cp_v::FT = constants.Cp_vapor #double check

    θ_dry,ρ_dry,T,p = calc_thermo_vars(θ_z,q_vap,ρ)

    r_vap::FT = q_vap / (1 - q_vap)
    R_m = R_v / (1 / r_vap + 1) + R_d / (1 + r_vap)
    cp_m = cp_v / (1 / r_vap + 1) + cp_d / (1 + r_vap)

    ρ_m = p / R_m / T

    return g / T * ρ_m * (R_m / cp_m - 1) / R_m
end


function dP_dz(P_z, ode_settings, z)

    FT = eltype(P_z)
    constants = ode_settings

    θ_z::FT = θl(z)
    q_vap::FT = qt(z)

    g::FT = constants.gconst
    R_d::FT = constants.Rd
    R_v::FT = constants.Rv
    cp_d::FT = constants.Cp_air
    cp_v::FT = constants.Cp_vapor #double check

    T::FT = in_situ_temperature(θ_z,q_vap,P_z,constants)
    T_virtual = T * (1 + 0.61 * q_vap)
    ρ::FT = P_z / (R_d * T_virtual)

    return -g * ρ
end


function calc_thermo_vars(θ_z,q_vap,ρ)
    θ_dry::FT = calc_θ_dry(constants,θ_z, q_vap)
    ρ_dry::FT = calc_ρ_dry_from_ρ(ρ, q_vap)
    T::FT = calc_T(constants,θ_dry, ρ_dry)
    p::FT = calc_p(constants,ρ_dry, T, q_vap)

    return θ_dry,ρ_dry,T,p
end
