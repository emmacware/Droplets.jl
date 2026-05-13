#check all functions for halo logic
#mpdata boundary condition?
#check top cell doesn't change, add MPDATA boundary condition

# condensation rate, radiation term condensation rate
# collision rates, radiative cooling rates,
# advection, diffusion,

function single_column_timestep(grid::scm_eulerian_arrays{FT}, dt::FT, droplets::droplet_attributes{FT},
        coagsettings::coag_settings{FT}, spatialsettings::spatial_settings{FT},
        condensationsettings::condensation_settings{FT}, coagdata::coagulation_run_spatial,
        conddata::condensation_data{FT}, raddata::radiation_data{FT},turbdata::turbulence_data{FT},
        diagnosticsettings::diagnostic_settings{FT},
        prescribed_w::Function, mpdatatmp::mpdata_tmp_1d, mpdatasettings::mpdata_settings_1d,
        constants::Constants{FT},
        scmsettings::scm_settings{FT},
        tkesettings::tke_settings{FT},
        i::Int

    ) where {FT<:AbstractFloat}
    output_idx = Int(floor(i*dt/spatialsettings.dt_output) +1)
    n_output = spatialsettings.dt_output / dt
    

    sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)
    # or i = 1
    if i % n_output == 0 || i == 1
        scm_fill_diagnostic_output(grid,coagdata,conddata,raddata,spatialsettings, output_idx)
    end

    if scmsettings.radiation_on
        radiation_function!(grid, spatialsettings, diagnosticsettings, constants, raddata, dt, i)
    end

    if scmsettings.condensation_on
        for t_cond in 1:scmsettings.n_cond
            grid.states.T_tmp .= T_from_theta.(grid.states.θ, grid.states.P, constants)
            condensation_time_step_spatial!(droplets, grid.states, nz, dt/scmsettings.n_cond, conddata, constants, condensationsettings, spatialsettings, raddata, i*dt)
        end
    end

    if scmsettings.coalescence_on && i*dt > scmsettings.spinup_time
        for t_coag in 1:scmsettings.n_coag
            coalescence_timestep!(scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)
        end
    end

    # Environmental advection
    mpdata_scm!(grid, dt, mpdatatmp, mpdatasettings, constants)

    update_droplet_positions!(droplets, prescribed_w, dt, spatialsettings, scmsettings, i)

    #change sorting indexes if droplets moved
    # droplets.I .= sort(droplets.I, by = k -> droplets.z_loc[k])
    sort!(droplets.I, by = k -> droplets.z_loc[k])

    recycle_precipitation!(droplets, grid, spatialsettings, diagnosticsettings,coagdata, constants, output_idx)
    # droplets.I .= sort(droplets.I, by = k -> droplets.cell_id[k])
    sort!(droplets.I, by = k -> droplets.cell_id[k])

    for g in eachindex(droplets.grid_range)
        start = findfirst(k -> droplets.cell_id[k] == g, coagdata.I)
        if start == nothing
            droplets.grid_range[g] = nothing
            continue
        end
        droplets.grid_range[g] = start:findlast(k -> droplets.cell_id[k] == g, droplets.I)
    end

    if scmsettings.turbulence_on
        turb_timestep!(grid, tkesettings, constants, dt, scmsettings, turbdata)
    end

    

    return nothing
end




function condensation_time_step_spatial!(droplets::droplet_attributes{FT}, state::scm_states{FT}, nz::Int, Δtg::FT, conddata::condensation_data{FT}, constants, 
    condsettings::condensation_settings{FT}, spatialsettings::spatial_settings{FT}, raddata::radiation_data{FT}, m_t::FT) where {FT<:AbstractFloat}
    # Calculate condensation time step for each droplet
    Ns = spatialsettings.Nz
    # FT = eltype(droplets.X)

    conddata.vol_change_helper .= 0
    # for idrop in eachindex(droplets.X)
    #     fill_cond_rad_term!(droplets,idrop,raddata,constants)
    # end

    Threads.@threads for z in 1:nz
    # for z in 1:nz
        if droplets.grid_range[z] == nothing
            continue
        end
        @views begin
        k = droplets.I[droplets.grid_range[z]]

        T_k = state.T_tmp[z]
        qv_k = state.qv[z]
        P_k = state.P[z]
        conddata.vol_change_helper[z] -= sum(droplets.X[k].* droplets.ξ[k])

        S_env = sat(qv_k, P_k) / esat(T_k)
        # if m_t<3600
        #     S_env = min(S_env, 1.01) # spinup supersat cap
        #     droplets.X[k] .= dXkappakohler_newtonraphson.(volume_to_radius.(droplets.X[k]), droplets.dry_r3[k], condsettings.kappa, T_k, S_env,constants, Δtg,10)#raddata.rad_term[k]) #giving up on DifferentialEquations for now, probably will put in more accurate integrator
        # else
        #     for idrop in k
        #         dXkappakohler_newtonraphson(droplets,idrop, condsettings.kappa, T_k, S_env,constants,raddata,absliq_r_interp, Δtg,10)#raddata.rad_term[k]) #giving up on DifferentialEquations for now, probably will put in more accurate integrator
        #     end   
        # end

        if m_t<3600
            S_env = min(S_env, 1.01) # spinup supersat cap
        end
        droplets.X[k] .= dXkappakohler_newtonraphson.(volume_to_radius.(droplets.X[k]), droplets.dry_r3[k], condsettings.kappa, T_k, S_env,constants, Δtg,10)#raddata.rad_term[k]) #giving up on DifferentialEquations for now, probably will put in more accurate integrator


        # dR_rad = drrad_term.(R, T_k, raddata.cond_rad_term[k], Δtg)
        # ratio = 0.0#dX_rad ./ (dX .+ dX_rad) # if they are opposite maybe need to handle differently? for now just a start, need to clean this anyways
        conddata.vol_change_helper[z] += sum(droplets.X[k].* droplets.ξ[k])
        dqv = conddata.vol_change_helper[z] .* constants.ρl ./ (state.ρ[z] .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
        conddata.condensation_src[z] += dqv
        state.qv[z] -= dqv
        state.T_tmp[z] += dqv * constants.L ./ (constants.Cp_air)
        state.θ[z] = theta_from_T(state.T_tmp[z], state.P[z], constants)
        state.ρ[z] = ρ_ideal_gas(state.P[z], state.T_tmp[z], state.qv[z], constants)
        dX_rad = raddata.cond_rad_term[k]
        # conddata.condensation_src[z] += vol_change_helper[z] .* constants.ρl ./ (state.ρ[z] .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)

        conddata.condensation_rad_net[z] += sum(dX_rad .* droplets.ξ[k]) .* constants.ρl ./ (state.ρ[z] .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
        conddata.condensation_rad_abs[z] += sum(abs.(dX_rad) .* droplets.ξ[k]) .* constants.ρl ./ (state.ρ[z] .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
        end # @views
    end
    
    # dqvap = - vol_change .* constants.ρl ./ (state.ρ .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
    # state.qv .+= dqvap
    # state.T_tmp .-= dqvap * constants.L ./ (constants.Cp_air)
    # state.θ .= theta_from_T(state.T_tmp, state.P, constants)
    # state.ρ .= ρ_ideal_gas(state.P, state.T_tmp, state.qv, constants)

    return
end



#currently handled in the turbulent vertical diffusion as surface fluxes
function update_surface_forcings!(grid, constants, scmsettings)
    ρ = ρ_calc_θ(grid.states.P[1], grid.states.θ[1], grid.states.qv[1], constants)
    Π = (grid.states.P[1] / constants.P0)^(constants.Rd / constants.Cp_air)

    grid.states.θ[1] += scmsettings.surface_sensible_heat_flux * scmsettings.Δt /
                         (ρ * constants.Cp_air * Π * grid.dz)

    mass_evap = scmsettings.surface_latent_heat_flux * scmsettings.Δt / (constants.L * grid.dz)
    grid.states.qv[1] += mass_evap / ρ

    grid.states.ρ[1] = ρ_calc_θ(grid.states.P[1], grid.states.θ[1], grid.states.qv[1], constants)
    return
end

function recycle_precipitation!(droplets, grid, spatialsettings, diagnosticsettings, coagdata, constants,output_step)
    # Simple precipitation recycling: if droplets are in the bottom cell and larger than a threshold, reinitialize them as aerosols in the top cell, adding mass to surface precipitation diagnostic
    start = findfirst(k -> droplets.z_loc[k] <= 0, coagdata.I) 
    if start == nothing
        return
    end
    @views begin
        fallendrops_idx = start:findlast(k -> droplets.z_loc[k] <= 0, coagdata.I)
        precip_mass = sum(droplets.X[fallendrops_idx] .* droplets.ξ[fallendrops_idx]) * constants.ρl
        grid.output.surface_precipitation[output_step] += precip_mass / (spatialsettings.area_per_grid * spatialsettings.dt_output) #
    end
    # Reinitialize precipitating droplets as aerosols in the top cell
    droplets.X[fallendrops_idx] .= 4pi/3 .*droplets.dry_r3[fallendrops_idx] # reset to dry mass
    droplets.cell_id[fallendrops_idx] .= spatialsettings.Nz # move to top
    droplets.z_loc[fallendrops_idx] .= spatialsettings.Z_max # move to top, should we move to below inversion?
    droplets.w_prime[fallendrops_idx] .= 0

    # droplets.ξ[fallendrops_idx] how to reset multiplicity? these are lots of initial drops that collided
    return
end
        

