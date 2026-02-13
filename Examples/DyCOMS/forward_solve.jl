#check all functions for halo logic
#mpdata boundary condition?
#check top cell doesn't change, add MPDATA boundary condition

function single_column_timestep(grid::scm_eulerian_arrays{FT}, dt::FT, droplets::droplet_attributes{FT},
        coagsettings::coag_settings{FT}, spatialsettings::spatial_settings{FT},
        condensationsettings::condensation_settings{FT},condensation_integrator, coagdata::coagulation_run_spatial, 
        diagnosticsettings::diagnostic_settings{FT},
        prescribed_w::Function, mpdata_tmp::mpdata_tmp_1d, mpdatasettings::mpdata_settings_1d,
        constants,
        scmsettings::scm_settings{FT}

    ) where {FT<:AbstractFloat}
    

    sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)

    #Update microphysics (condensation, coagulation)
    condensation_time_step_spatial!(droplets, grid.states,nz, dt, condensation_integrator, constants,condensationsettings)
    coalescence_timestep!(scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)

    #Droplet motion (advection and settling)
    update_droplet_positions!(droplets, prescribed_w, dt, spatialsettings)
    
    # Environmental advection 
    mpdata_scm!(grid, dt, mpdata_tmp, mpdatasettings,constants)

    #surface forcings
    update_surface_forcings!(grid, constants, scmsettings)

    return nothing
end


function condensation_time_step_spatial!(droplets, state,nz, Δtg,condensation,constants,condsettings)
    # Calculate the condensation time step for each droplet
    droplets.X
    lnR = log.(volume_to_radius.(droplets.X))

    condensation.u.lnR .= lnR
    condensation.u.qvap .= state.qv
    condensation.u.T .= state.T

    condensation.p = (nz,droplets, state, constants, condsettings)
    
    step!(condensation, Δtg, true)

    droplets.X .= radius_to_volume.(exp.(condensation.u.lnR))
    state.qv .= condensation.u.qvap
    state.T .= condensation.u.T
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
