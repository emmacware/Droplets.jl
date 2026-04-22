#check all functions for halo logic
#mpdata boundary condition?
#check top cell doesn't change, add MPDATA boundary condition

# condensation rate, radiation term condensation rate
# collision rates, radiative cooling rates,
# advection, diffusion,

function single_column_timestep(grid::scm_eulerian_arrays{FT}, dt::FT, droplets::droplet_attributes{FT},
        coagsettings::coag_settings{FT}, spatialsettings::spatial_settings{FT},
        condensationsettings::condensation_settings{FT}, coagdata::coagulation_run_spatial,
        conddata::condensation_data{FT}, raddata::radiation_data{FT},
        diagnosticsettings::diagnostic_settings{FT},
        prescribed_w::Function, mpdata_tmp::mpdata_tmp_1d, mpdatasettings::mpdata_settings_1d,
        constants,
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

    radiation_function!(grid,spatialsettings, diagnosticsettings, constants,raddata, dt, i)

    for t_cond in 1:10
        # println("Condensation substep: ", t_cond)
        grid.states.T_tmp .= T_from_theta.(grid.states.θ,grid.states.P,constants)
        condensation_time_step_spatial!(droplets, grid.states,nz, dt/10, conddata, constants,condensationsettings,spatialsettings,coagdata,raddata,i)
    end
    if i > 3600
        coalescence_timestep!(scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)
    end
    #Droplet motion (advection and settling)
    update_droplet_positions!(droplets, prescribed_w, dt, spatialsettings)
    # Environmental advection 
    mpdata_scm!(grid, dt, mpdata_tmp, mpdatasettings,constants)

    # duynkerke_tke_timestep!(grid, tkesettings, constants, dt)
    no_prog_tke_timestep!(grid, tkesettings, constants, dt)
    # smag_lilly_timestep!(grid, tkesettings, constants, dt)
    # smag_lilly_ca_timestep!(grid, tkesettings, constants, dt)
    turbulent_droplet_diffusion!(droplets, grid, tkesettings, dt,coagdata,constants)
    
    #change sorting indexes if droplets moved
    coagdata.I .= sort(coagdata.I, by = k -> droplets.z_loc[k])
    recycle_precipitation!(droplets, grid, spatialsettings, diagnosticsettings, coagdata, constants,output_idx)
    coagdata.I .= sort(coagdata.I, by = k -> droplets.cell_id[k])
    for g in eachindex(droplets.grid_range)
        start = findfirst(k -> droplets.cell_id[k] == g, coagdata.I) 
        if start == nothing
            droplets.grid_range[g] = nothing
            continue
        end
        droplets.grid_range[g] = start:findlast(k -> droplets.cell_id[k] == g, coagdata.I)
    end

    
    #surface forcings
    # update_surface_forcings!(grid, constants, scmsettings)


    return nothing
end



function condensation_time_step_spatial!(droplets, state, nz, Δtg,condensation,constants,condsettings,spatialsettings,coagdata,raddata,m_t)
    # Calculate condensation time step for each droplet
    Ns = length(droplets.X)
    FT = eltype(droplets.X)

    R_array = range(lookup_lw_cld.bounds[1],lookup_lw_cld.bounds[2],lookup_lw_cld.dims[3])
    vol_change = zeros(FT, nz)

    for z in 1:nz
        if droplets.grid_range[z] == nothing
            continue
        end
        k = coagdata.I[droplets.grid_range[z]]

        X = droplets.X[k]
        R = volume_to_radius.(droplets.X[k])
        T_k = state.T_tmp[z]
        qv_k = state.qv[z]
        P_k = state.P[z]

        # # rad_term = zeros(size(X))
        # for idrop in 1:length(R)
        #     if R[idrop] >= 2.5e-6 && X[idrop] >= diagnosticsettings.aerosol_cloud_cuttoff && X[idrop] < diagnosticsettings.cloud_rain_cuttoff
        #         for ibnd in 1:16
        #             fnd_k = raddata.flux_net_droplet[z,ibnd]
        #             extliq,ssaliq,asyliq = LookUpTables.getview_liqdata(lookup_lw_cld,ibnd)
        #             absliq = extliq.*(1.0.-ssaliq)                   

        #             abs_drop = RRTMGP.Optics.interp1d_equispaced(R[idrop]*1.0e6,R_array,absliq)
        #             Qa = abs_drop * X[idrop] * constants.ρl * 1000.0 / (pi *R[idrop]^2)

        #             #if R[idrop]>10.0e-6 
        #             #    println(R[idrop]*1.0e6,' ',Qa)
        #             #end
        #             raddata.cond_rad_term[k[idrop]] += Qa * raddata.flux_net_droplet[z,ibnd]
        #         end
        #     else
        #         raddata.cond_rad_term[k[idrop]] = 0.0
        #     end
        # end

        S_env = sat(qv_k, P_k) / esat(T_k)
        if m_t<2400
            S_env = min(S_env, 1.01) # spinup supersat cap
        end
        # dR = drkappakohler.(R, droplets.dry_r3[k], condsettings.kappa, T_k, S_env, Δtg,rad_term=rad_term[k]) #giving up on DifferentialEquations for now, probably will put in more accurate integrator
        # dR = drkappakohler.(R, droplets.dry_r3[k], condsettings.kappa, T_k, S_env, Δtg,rad_term=0.0)#raddata.rad_term[k]) #giving up on DifferentialEquations for now, probably will put in more accurate integrator
        dX = dXkappakohler_newtonraphson.(R, droplets.dry_r3[k], condsettings.kappa, T_k, S_env, Δtg,10)#raddata.rad_term[k]) #giving up on DifferentialEquations for now, probably will put in more accurate integrator
        # dR_rad = drrad_term.(R, T_k, raddata.cond_rad_term[k], Δtg)
        # dX = (4 * pi .* R.^2 .* dR) .* Δtg
        # dX_rad = (4 * pi .* R.^2 .* dR_rad) .* Δtg
        dry_vol = radius_to_volume.((droplets.dry_r3[k]).^(1/3))
        ratio = 0.0#dX_rad ./ (dX .+ dX_rad) # if they are opposite maybe need to handle differently? for now just a start, need to clean this anyways
        dX .= max.(dX, dry_vol .- droplets.X[k]) # snap to dry mass if smaller
        # dX_rad .= dX .* ratio # calculate how much of the total growth is due to radiation, for diagnostics
        vol_change[z] += sum(dX .* droplets.ξ[k])
        conddata.condensation_src[z] += vol_change[z] .* constants.ρl ./ (state.ρ[z] .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
        # conddata.condensation_rad_net[z] += sum(dX_rad .* droplets.ξ[k]) .* constants.ρl ./ (state.ρ[z] .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
        # conddata.condensation_rad_abs[z] += sum(abs.(dX_rad) .* droplets.ξ[k]) .* constants.ρl ./ (state.ρ[z] .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
        droplets.X[k] .+= dX

    end
    
    dqvap = - vol_change .* constants.ρl ./ (state.ρ .* spatialsettings.area_per_grid * spatialsettings.z_grid_height)
    state.qv .+= dqvap
    state.T_tmp .-= dqvap * constants.L ./ (constants.Cp_air)
    state.θ .= theta_from_T(state.T_tmp, state.P, constants)
    state.ρ .= ρ_ideal_gas(state.P, state.T_tmp, state.qv, constants)

    return
end


# function update_surface_forcings!(grid, constants, scmsettings)
#     grid.states.ρ[1] = ρ_calc_θ(grid.states.P[1],grid.states.θ[1],grid.states.qv[1],constants)
#     dT = (scmsettings.surface_sensible_heat_flux * scmsettings.Δt / (grid.states.ρ[1] * constants.Cp_air * grid.dz))
#     dθ = dT.* (constants.P0 ./ grid.states.P[1]).^(constants.Rd / constants.Cp_air)
#     grid.states.θ[1] += dθ #theta_from_T(grid.states.T_tmp[1], grid.states.P[1], constants)
#     # grid.states.ρ[1] = ρ_calc_θ(grid.states.P[1],grid.states.θ[1],grid.states.qv[1],constants)
#     mass_evap = scmsettings.surface_latent_heat_flux * scmsettings.Δt / (constants.L * grid.dz)
#     mass_per_volume = grid.states.ρ[1]
#     mass_vapor = grid.states.qv[1] * mass_per_volume
#     grid.states.qv[1] = (mass_vapor + mass_evap) / (mass_per_volume + mass_evap) # update qv based on added vapor mass, assuming it displaces an equal mass of air (so total mass in cell increases by mass_evap)    #update ρ and theta after surface fluxes
#     grid.states.ρ[1] = ρ_calc_θ(grid.states.P[1],grid.states.θ[1],grid.states.qv[1],constants)
#     return
# end

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
    fallendrops_idx = start:findlast(k -> droplets.z_loc[k] <= 0, coagdata.I)
    precip_mass = sum(droplets.X[fallendrops_idx] .* droplets.ξ[fallendrops_idx]) * constants.ρl
    grid.output.surface_precipitation[output_step] += precip_mass / (spatialsettings.area_per_grid * spatialsettings.dt_output) #

    # Reinitialize precipitating droplets as aerosols in the top cell
    droplets.X[fallendrops_idx] .= droplets.dry_r3[fallendrops_idx] # reset to dry mass
    droplets.cell_id[fallendrops_idx] .= spatialsettings.Nz # move to top
    droplets.z_loc[fallendrops_idx] .= spatialsettings.Z_max # move to top, should we move to below inversion?
    # droplets.ξ[fallendrops_idx] how to reset multiplicity? these are lots of initial drops that collided
    return
end
        

