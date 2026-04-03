#check all functions for halo logic
#mpdata boundary condition?
#check top cell doesn't change, add MPDATA boundary condition

function single_column_timestep(grid::scm_eulerian_arrays{FT}, dt::FT, droplets::droplet_attributes{FT},
        coagsettings::coag_settings{FT}, spatialsettings::spatial_settings{FT},
        condensationsettings::condensation_settings{FT},condensation_integrator, coagdata::coagulation_run_spatial, 
        diagnosticsettings::diagnostic_settings{FT},
        prescribed_w::Function, mpdata_tmp::mpdata_tmp_1d, mpdatasettings::mpdata_settings_1d,
        constants,
        scmsettings::scm_settings{FT},
        tkesettings::tke_settings{FT},

    ) where {FT<:AbstractFloat}
    
    

    sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)

    println("radiation")
    flux_net_droplet = radiation_function!(grid,spatialsettings, diagnosticsettings, constants, dt)

    #Update microphysics (condensation, coagulation)
    for t_cond in 1:100
        # println("Condensation substep: ", t_cond)
        condensation_time_step_spatial!(droplets, grid.states,nz, dt/100, condensation_integrator, constants,condensationsettings,spatialsettings)
    end
    coalescence_timestep!(scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)
    #Droplet motion (advection and settling)
    update_droplet_positions!(droplets, prescribed_w, dt, spatialsettings)
    
    # Environmental advection 
    mpdata_scm!(grid, dt, mpdata_tmp, mpdatasettings,constants)

    #surface forcings
    update_surface_forcings!(grid, constants, scmsettings)


    return nothing
end


function condensation_time_step_spatial!(droplets, state,nz, Δtg,condensation,constants,condsettings,spatialsettings)
    # Calculate the condensation time step for each droplet
    # droplets.X
    # lnr = log.((volume_to_radius.(droplets.X)))
    Ns = length(droplets.X)
    FT = eltype(droplets.X)

    
    mass_change = zeros(FT, nz)

    for k in 1:Ns
        z = droplets.cell_id[k]
        if z < 1
            continue
        end
        
        R = volume_to_radius(droplets.X[k])
        T_k = state.T[z]
        qv_k = state.qv[z]
        P_k = state.P[z]

        S_env = sat(qv_k, P_k) / esat(T_k)
        dR = drkappakohler(R, droplets.dry_r3[k], condsettings.kappa, T_k, S_env, Δtg)
        dX = (4 * pi * R^2 * dR) * Δtg
        if droplets.X[k] + dX < 0
            dX = -droplets.X[k] + radius_to_volume((droplets.dry_r3[k])^(1/3))  # prevent negative mass
        end

        mass_change[z] += dX * droplets.ξ[k]
        droplets.X[k] += dX



        # mass_change[z] += (radius_to_volume(exp(final_lnr)) - droplets.X[k]) * droplets.ξ[k]
        # droplets.X[k] = radius_to_volume(exp(final_lnr))
    end
    
    dqvap = - mass_change .* constants.ρl ./ (state.ρ .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
    state.qv .+= dqvap
    state.T .-= dqvap * constants.L ./ (constants.Cp_air)
    state.θ .= theta_from_T(state.T, state.P, constants)
    state.ρ .= ρ_ideal_gas(state.P, state.T, state.qv, constants)

    return
end



function update_surface_forcings!(grid, constants, scmsettings)
    grid.states.T[1] += scmsettings.surface_sensible_heat_flux * scmsettings.Δt / (grid.states.ρ[1] * constants.Cp_air * grid.dz)
    grid.states.θ[1] = theta_from_T(grid.states.T[1], grid.states.P[1], constants)
    grid.states.ρ[1] = ρ_ideal_gas(grid.states.P[1], grid.states.T[1], grid.states.qv[1], constants)
    grid.states.qv[1] += scmsettings.surface_latent_heat_flux * scmsettings.Δt / (grid.states.ρ[1] * constants.L * grid.dz)
    return
end
