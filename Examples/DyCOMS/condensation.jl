using ComponentArrays

function condensation_time_step_spatial(droplets, state,nz, Δtg,condensation,constants)
    # Calculate the condensation time step for each droplet
    droplets.X
    lnR = log.(volume_to_radius.(droplets.X))

    condensation.u.lnR .= lnR
    condensation.u.qvap .= state.qv
    condensation.u.T .= state.T

    condensation.p = (nz,droplets, state, constants)
    
    step!(condensation, Δtg, true)

    droplets.X .= radius_to_volume.(exp.(condensation.u.lnR))
    state.qv .= condensation.u.qvap
    state.T .= condensation.u.T
    state.θ .= state.T .* (1000 ./ state.P).^(constants.Rd / constants.Cp_air)
    state.ρ .= state.P ./ (constants.Rd .* state.T .* (1 .+ 0.61 .* state.qv))

    return
end


function condensation_rhs(du, u, p, t)
    du .= zero.(u)  # Initialize du to zero, preserving the structure and types of u
    FT = eltype(u)
    lnR = u.lnR
    qvap_col = u.qvap
    T_col = u.T
    nz,drops,grid,constants = p
    R = exp.(lnR)
    T_v = T_col .* (1 .+ 0.61 .* qvap_col)
    Rd = constants.Rd

    for z in 1:length(nz)
        k = z + 1 #account for halo cell
        P::FT = grid.P[k]
        T::FT  = T_col[k]
        qv::FT = qvap_col[k]
        Tv::FT = T_v[k]
        ρ_air::FT = P / Rd / Tv
        S_env::FT = sat.(qv, P) ./ esat(T)
        #find all indexes where drops.cell_id == k
        R_idx = findall(drops.cell_id .== k)

        dR = drkappakohler.(R[R_idx],drops.dry_r3[R_idx],drops.κ[R_idx],T,S_env,t)
        dX = 4 * π .* R[R_idx].^2 .* dR

        dqvap = - sum(dX .* drops.ξ[R_idx] .* constants.ρl / ρ_air)
        dT = - dqvap * constants.L / constants.Cp_air

        du.lnR[R_idx] .= FT.(dR ./ R[R_idx])
        du.qvap[k] = FT(dqvap)
        du.T[k] = FT(dT)
    end
end

# function condensation_rhs_single_cell(du, u, p, t)
#     lnR = u.lnR
#     du.lnR .= 0
#     qv = u.qvap
#     T = u.T
#     drops,grid,constants,k = p
#     R = exp.(lnR)
#     Tv = T .* (1 + 0.61 * qvap)
#     Rd = constants.Rd

#     P = grid.P[k]
#     ρ_air = P / Rd / Tv
#     S_env = sat.(qv, P) ./ esat(T)

#     R_idx = findall(drops.cell_id .== k)

#     dR = drkappakohler(R[R_idx],drops.dry_r3[R_idx],drops.kappa[R_idx],T,Senv,t)
#     dX = 4 * π .* R[R_idx]^2 * dR

#     dqvap = - sum(dX .* drops.ξ[R_idx] .* ρₗ / ρ_air)
#     dT = - dqvap * constants.L / constants.Cp_air

#     du.lnR[R_idx] = dR ./ R[R_idx]
#     du.qvap = dqvap
#     du.T = dT

# end