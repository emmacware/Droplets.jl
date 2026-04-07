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
        i::Int

    ) where {FT<:AbstractFloat}
    
    

    sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)

    println("radiation")
    flux_net_droplet = radiation_function!(grid,spatialsettings, diagnosticsettings, constants, dt, i)

    #Update microphysics (condensation, coagulation)
    for t_cond in 1:100
        # println("Condensation substep: ", t_cond)
        condensation_time_step_spatial!(droplets, grid.states,nz, dt/100, condensation_integrator, constants,condensationsettings,spatialsettings,coagdata,flux_net_droplet)
    end
    coalescence_timestep!(scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)
    #Droplet motion (advection and settling)
    update_droplet_positions!(droplets, prescribed_w, dt, spatialsettings)

    # Environmental advection 
    mpdata_scm!(grid, dt, mpdata_tmp, mpdatasettings,constants)
    tke_timestep!(grid, tkesettings, constants, dt)
    turbulent_droplet_diffusion!(droplets, grid, tkesettings, dt,coagdata)
    
    #change sorting indexes if droplets moved
    coagdata.I .= sort(coagdata.I, by = i -> droplets.cell_id[i])
    for g in eachindex(droplets.grid_range)
        start = findfirst(i -> droplets.cell_id[i] == g, coagdata.I) 
        if start == nothing
            droplets.grid_range[g] = nothing
            continue
        end
        droplets.grid_range[g] = start:findlast(i -> droplets.cell_id[i] == g, coagdata.I)
    end
    
    #surface forcings
    update_surface_forcings!(grid, constants, scmsettings)


    return nothing
end



function condensation_time_step_spatial!(droplets, state,nz, Δtg,condensation,constants,condsettings,spatialsettings,coagdata,flux_net_droplets)
    # Calculate condensation time step for each droplet
    Ns = length(droplets.X)
    FT = eltype(droplets.X)

    
    mass_change = zeros(FT, nz)

    for z in 1:nz
        if droplets.grid_range[z] == nothing
            continue
        end
        k = coagdata.I[droplets.grid_range[z]]

        R = volume_to_radius.(droplets.X[k])
        T_k = state.T[z]
        qv_k = state.qv[z]
        P_k = state.P[z]
        #fnd_k = flux_net_droplets[k,:]

        S_env = sat(qv_k, P_k) / esat(T_k)
        dR = drkappakohler.(R, droplets.dry_r3[k], condsettings.kappa, T_k, S_env, Δtg) #giving up on DifferentialEquations for now, probably will put in more accurate integrator
        dX = (4 * pi .* R.^2 .* dR) .* Δtg
        dry_vol = radius_to_volume.((droplets.dry_r3[k]).^(1/3))
        dX .= max.(dX, dry_vol .- droplets.X[k]) # snap to dry mass if smaller
        mass_change[z] += sum(dX .* droplets.ξ[k])
        droplets.X[k] .+= dX

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
    #update ρ and theta after surface fluxes
    grid.states.θ[1] = theta_from_T(grid.states.T[1], grid.states.P[1], constants)
    grid.states.ρ[1] = ρ_ideal_gas(grid.states.P[1], grid.states.T[1], grid.states.qv[1], constants)
    return
end