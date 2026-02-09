

function single_column_timestep(grid::scm_eulerian_arrays{FT}, dt::FT, droplets::droplet_attributes{FT},
        coagsettings::coag_settings{FT}, spatial::spatial_settings{FT},
        condensation::condensation_settings{FT},condensation_integrator, coag_data::coagulation_run_spatial, 
        diagnosticsettings::diagnostic_settings{FT},
        w_function::Function, tmp::mpdata_tmp_1d, mpdatasettings::mpdata_settings_1d,
        constants,
        scmsettings::scm_settings{FT}

    ) where {FT<:AbstractFloat}
    
    #Update thermodynamcis

    #Droplet radiation
    sd_fill_diagnostics(droplets, grid, spatial, diagnosticsettings)

    #Update microphysics (condensation, coagulation)
    condensation_time_step_spatial(droplets, grid.states,nz, dt, condensation_integrator, constants)
    coalescence_timestep!(scmsettings.coag_threading, scmsettings.scheme, droplets, coag_data, coagsettings)

    #Droplet sedimentation
    update_droplet_positions(droplets, w_function, dt)
    
    # Advection function
    mpdata_scm(grid, dt, tmp, mpdatasettings,constants)

    return nothing
end