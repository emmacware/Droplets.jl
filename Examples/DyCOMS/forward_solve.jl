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
        absliq_r_interp,
        i::Int

    ) where {FT<:AbstractFloat}
    output_idx = Int(floor(i*dt/spatialsettings.dt_output) +1)
    n_output = spatialsettings.dt_output / dt
    

    sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)
    # or i = 1
    if i % n_output == 0 || i == 1
        scm_fill_diagnostic_output(grid,coagdata,conddata,raddata,spatialsettings,constants, output_idx)
    end


    radiation_function!(scmsettings.radiation,grid, spatialsettings, diagnosticsettings, constants, raddata, dt, i)


    T_from_theta!(grid.states.T_tmp,grid.states.θ, grid.states.P, constants)
    for _ in 1:scmsettings.n_cond
        condensation_time_step_spatial!(scmsettings.condensation,droplets, grid.states, nz, dt/scmsettings.n_cond, conddata, constants, condensationsettings, spatialsettings, raddata,scmsettings,absliq_r_interp, i*dt)
    end



    for _ in 1:scmsettings.n_coag
        coalescence_timestep!(scmsettings.coalescence,scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)
    end


    # Environmental advection
    mpdata_scm!(scmsettings.advection,grid, dt, mpdatatmp, mpdatasettings, constants)
    turb_timestep!(scmsettings.turbulence,grid, tkesettings, constants, dt, scmsettings, turbdata)

    update_droplet_positions!(scmsettings.motion,scmsettings.advection,scmsettings.settling,droplets, prescribed_w, dt, spatialsettings, scmsettings, i)
    sort!(droplets.I, by = k -> droplets.z_loc[k])
    recycle_precipitation!(scmsettings.recycling,droplets, grid, spatialsettings, diagnosticsettings,coagdata, constants, output_idx)
    sort!(droplets.I, by = k -> droplets.cell_id[k])

    fill!(droplets.grid_range, 1:0)
    nsd = length(droplets.I)
    i = 1
    while i <= nsd
        g  = droplets.cell_id[droplets.I[i]]
        
        lo = i
        while i <= nsd && droplets.cell_id[droplets.I[i]] == g
            i += 1
        end
        g < 1 && continue
        droplets.grid_range[g] = lo:(i - 1)
    end

    # sort!(droplets.I, by = k -> droplets.cell_id[k])

    # for g in eachindex(droplets.grid_range)
    #     start = findfirst(k -> droplets.cell_id[k] == g, coagdata.I)
    #     if start == nothing
    #         droplets.grid_range[g] = 1:0
    #         continue
    #     end
    #     droplets.grid_range[g] = start:findlast(k -> droplets.cell_id[k] == g, droplets.I)
    # end


    return nothing
end


function recycle_precipitation!(::DynOFF,droplets, grid, spatialsettings, diagnosticsettings, coagdata, 
    constants,output_step)
end

function recycle_precipitation!(::DynON, droplets, grid, spatialsettings, diagnosticsettings, coagdata, constants, output_step)
    FT = eltype(droplets.X)
    precip_mass = zero(FT)
    for k in droplets.I
        droplets.z_loc[k] > 0 && continue
        # droplets.z_loc[k] == -10 && continue # already recycled

        #dont let turbulence cause "precipitation"
        if droplets.X[k] < diagnosticsettings.aerosol_cloud_cuttoff
            droplets.z_loc[k] = rand() * spatialsettings.z_grid_height
            droplets.cell_id[k] = 1
            continue
        end

        precip_mass += droplets.X[k] * droplets.ξ[k] * constants.ρl
        # droplets.X[k] = 0#FT(4π/3) * droplets.dry_r3[k]
        eligible = findall(r -> length(r) >= 2, droplets.grid_range[1:end-1])
        recycle_idx = isempty(eligible) ? argmax(length.(droplets.grid_range[1:end-1])) : eligible[argmin(length.(droplets.grid_range[eligible]))]
        drop_to_split = droplets.I[droplets.grid_range[recycle_idx][end]]

        droplets.ξ[k] = floor(Int,droplets.ξ[drop_to_split]/2)
        droplets.X[k] = droplets.X[drop_to_split]
        droplets.dry_r3[k] = droplets.dry_r3[drop_to_split]
        droplets.cell_id[k] = droplets.cell_id[drop_to_split]
        droplets.z_loc[k] = droplets.z_loc[drop_to_split]
        droplets.w_prime[k] = FT(0.0)
        droplets.ξ[drop_to_split] -= droplets.ξ[k]
    end
    grid.output.surface_precipitation[output_step] += precip_mass  / (spatialsettings.area_per_grid * spatialsettings.dt_output)
    return
end
        

# function recycle_precipitation!(::DynON,droplets, grid, spatialsettings, diagnosticsettings, coagdata, constants,output_step)
#     # Simple precipitation recycling: if droplets are in the bottom cell and larger than a threshold, reinitialize them as aerosols in the top cell, adding mass to surface precipitation diagnostic
#     start = findfirst(k -> droplets.z_loc[k] <= 0, coagdata.I) 
#     if start == nothing
#         return
#     end
#     fallendrops_idx = start:findlast(k -> droplets.z_loc[k] <= 0, coagdata.I)
#     precip_mass = sum(droplets.X[fallendrops_idx] .* droplets.ξ[fallendrops_idx]) * constants.ρl
#     grid.output.surface_precipitation[output_step] += precip_mass / (spatialsettings.area_per_grid * spatialsettings.dt_output) #
#     # Reinitialize precipitating droplets as aerosols in the top cell
#     # droplets.X[fallendrops_idx] .= 4pi/3 .*droplets.dry_r3[fallendrops_idx] # reset to dry mass
#     # droplets.cell_id[fallendrops_idx] .= spatialsettings.Nz # move to top
#     # droplets.z_loc[fallendrops_idx] .= spatialsettings.Z_max # move to top, should we move to below inversion?
#     # droplets.w_prime[fallendrops_idx] .= 0
#     droplets.ξ[fallendrops_idx] .= 0 # set multiplicity to 0 so they don't affect diagnostics until they grow again, but keep mass so they can
#     # split_highest_multiplicity!(droplets) # split the highest multiplicity drops into many smaller drops to avoid having a few huge drops with large multiplicity that dominate the diagnostics. This is a bit hacky but it allows us to keep the mass of the precipitating drops without having a few huge drops with large multiplicity that dominate the diagnostics. We could also just set the mass to 0 and let them grow again, but this way we can keep track of the precipitating mass in the diagnostics without having a few huge drops with large multiplicity that dominate the diagnostics.

#     # droplets.ξ[fallendrops_idx] how to reset multiplicity? these are lots of initial drops that collided
#     return
# end
        