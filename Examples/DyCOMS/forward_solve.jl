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


    T_from_theta!(grid.states.T_tmp,grid.states.θ, grid.states.P, grid.states.qv,constants)
    for _ in 1:scmsettings.n_cond
        condensation_time_step_spatial!(scmsettings.condensation,droplets, grid.states, nz, dt/scmsettings.n_cond, conddata, constants, condensationsettings, spatialsettings, raddata,scmsettings,absliq_r_interp, i*dt)
    end



    for _ in 1:scmsettings.n_coag
        coalescence_timestep!(scmsettings.coalescence,scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings)
    end


    # Environmental advection
    turb_timestep!(scmsettings.turbulence,grid, tkesettings, constants, dt, scmsettings, turbdata)
    compute_ql_at_cell!.(grid.states, 1:nz,constants)
    qt = grid.states.qv .+ grid.states.ql_tmp
    θ_l = θl.(grid.states.P, grid.states.T_tmp, grid.states.ql_tmp, grid.states.qv, constants)
    fill_grid_ranges!(droplets)
    compute_ql_at_cell!.(grid.states, 1:nz,constants)
    grid.states.qv .= qt .- grid.states.ql_tmp
    grid.states.θ .= θ_from_θl.(grid.states.P, θ_l, grid.states.ql_tmp, grid.states.qv, constants)
    mpdata_scm!(scmsettings.advection, scmsettings.thermo_feedback, grid, dt, mpdatatmp, mpdatasettings, constants)
    update_droplet_positions!(scmsettings.motion,scmsettings.advection,scmsettings.settling,droplets, prescribed_w, dt, spatialsettings, scmsettings, i)

    sort!(droplets.I, by = k -> droplets.z_loc[k])
    recycle_precipitation!(scmsettings.recycling,droplets, grid, spatialsettings, diagnosticsettings,coagdata, constants, output_idx)
    recycle_top_escape!(scmsettings.top_escape, droplets, spatialsettings)
    sort!(droplets.I, by = k -> droplets.cell_id[k])
    # keep_layer_filled!(droplets, grid, spatialsettings, diagnosticsettings, constants)
    # sort!(droplets.I, by = k -> droplets.cell_id[k])

    fill_grid_ranges!(droplets)

    return nothing
end


function recycle_precipitation!(::DynOFF,droplets, grid, spatialsettings, diagnosticsettings, coagdata, 
    constants,output_step)
end

function recycle_precipitation!(::DynON, droplets, grid, spatialsettings, diagnosticsettings, coagdata, constants, output_step)
    FT = eltype(droplets.X)
    precip_mass = zero(FT)

    # grid_range must be consistent with cell_id order; the caller sorted by z_loc,
    # so refresh here before using grid_range to find drop_to_split
    fill_grid_ranges!(droplets)

    nz = length(droplets.grid_range)
    raw_inv    = findfirst(k -> grid.states.qv[k] <= FT(0.008), 1:nz-1)
    z_inv_idx  = raw_inv === nothing ? nz - 1 : raw_inv - 1
    raw_400    = findfirst(k -> grid.centers_z[k] >= 400, 1:nz-1)
    idx_400m   = raw_400 === nothing ? 1 : raw_400
    cloud_cells = filter(k -> !isempty(droplets.grid_range[k]), idx_400m:z_inv_idx)
    recycle_idx = isempty(cloud_cells) ? z_inv_idx :
                  argmin(k -> length(droplets.grid_range[k]), cloud_cells)

    for k in droplets.I
        droplets.z_loc[k] > 0 && droplets.ξ[k] > 0 && continue

        #dont let turbulence cause "precipitation"
        if droplets.X[k] < diagnosticsettings.aerosol_cloud_cuttoff
            droplets.z_loc[k] = rand() * spatialsettings.z_grid_height
            droplets.cell_id[k] = 1
            continue
        end

        precip_mass += droplets.X[k] * droplets.ξ[k] * constants.ρl
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

function recycle_top_escape!(::DynOFF, droplets, spatialsettings)
end

function recycle_top_escape!(::DynON, droplets, spatialsettings)
    FT = eltype(droplets.X)
    Z_max = spatialsettings.Z_max
    T = T_from_theta(grid.states.θ[1], grid.states.P[1],grid.states.qv[1], constants)
    qv_k = grid.states.qv[1]
    P_k = grid.states.P[1]
    S_env = sat(qv_k,P_k)/esat(T)
    if S_env > 0.95
        S_env = 0.95
    end

    for k in droplets.I
        droplets.z_loc[k] < Z_max && continue
        droplets.z_loc[k] = FT(0.1)
        droplets.cell_id[k] = 1
        find_equilibrium_radius.(droplets,k, 0.95, T, S_env,constants)
    end

    return
end

function fill_grid_ranges!(droplets)
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
        g > length(droplets.grid_range) && continue
        droplets.grid_range[g] = lo:(i - 1)
    end
end

function keep_layer_filled!(droplets, grid, spatialsettings, diagnosticsettings, constants; min_count::Int = 400)
    FT = eltype(droplets.X)
    fill_grid_ranges!(droplets)  # ensure ranges are current
    nz = length(droplets.grid_range)

    dz = spatialsettings.z_grid_height
    z_inv_idx = findfirst(k -> grid.states.qv[k] <= FT(0.008), 1:length(droplets.grid_range)-1) - 1
    idx_400m = findfirst(k -> grid.centers_z[k] >= 400, 1:length(droplets.grid_range)-1)

    all(c -> length(droplets.grid_range[c]) >= min_count, idx_400m:z_inv_idx) && return nothing
    cell_sd = [collect(droplets.I[droplets.grid_range[g]]) for g in 1:nz]


    for c in idx_400m:z_inv_idx
        N_c = length(cell_sd[c])
        N_c >= min_count && continue
        # println("Cell ", c, " has only ", N_c, " SDs; needs ", min_count)
        N_c == 0 && continue
        N_needed = min_count - N_c + 5

        # donor = cell with most SDs excluding c
        donor = 0
        donor_count = 0
        for g in 1:idx_400m
            g == c && continue
            lg = length(cell_sd[g])
            if lg > donor_count
                donor_count = lg
                donor = g
            end
        end

        # println("Donor cell: ", donor, " with count: ", donor_count)

        # Check if donor is valid
        if donor == 0 || donor_count <= 40
            println("No valid donor found or donor count too low.")
            # N_needed = 5
            break
        end

        # zero-ξ donor slots are already vacated — use them first at no cost
        d_idx_zero = filter(i -> droplets.ξ[i] == 0, cell_sd[donor])
        d_idx_active = filter(i -> droplets.ξ[i] > 0, cell_sd[donor])
        isempty(d_idx_active) && isempty(d_idx_zero) && continue

        n_free = length(d_idx_zero)
        n_merge_needed = max(0, N_needed - n_free)

        # sort active SDs by radius; merge tightest adjacent pairs to free more slots
        d_radii = volume_to_radius.(droplets.X[d_idx_active])
        order = sortperm(d_radii)
        d_idx_s = d_idx_active[order]
        n_active = length(d_idx_s)

        free_pos = Int[]
        keep_pos = Int[]
        if n_merge_needed > 0 && n_active >= 2
            diffs = diff(volume_to_radius.(droplets.X[d_idx_s]))
            gap_order = sortperm(diffs)
            used = falses(n_active)
            for g in gap_order
                length(free_pos) == n_merge_needed && break
                if !used[g] && !used[g+1]
                    push!(keep_pos, g)
                    push!(free_pos, g + 1)
                    used[g] = used[g+1] = true
                end
            end
        end

        # exact ξ conservation: absorb freed SD's multiplicity into its pair partner
        for (kp, fp) in zip(keep_pos, free_pos)
            ξ_sum = droplets.ξ[d_idx_s[kp]] + droplets.ξ[d_idx_s[fp]]
            droplets.X[d_idx_s[kp]] = (droplets.X[d_idx_s[kp]] * droplets.ξ[d_idx_s[kp]] + droplets.X[d_idx_s[fp]] * droplets.ξ[d_idx_s[fp]]) / ξ_sum
            droplets.ξ[d_idx_s[kp]] = ξ_sum
            droplets.ξ[d_idx_s[fp]] = 0
        end

        merged_free = d_idx_s[free_pos]
        free_set = Set(free_pos)
        n_zero_borrowed = min(N_needed, n_free)
        borrowed = vcat(d_idx_zero[1:n_zero_borrowed], merged_free)
        remaining_active = d_idx_s[[i for i in 1:n_active if i ∉ free_set]]
        remaining = vcat(remaining_active, d_idx_zero[n_zero_borrowed+1:end])

        # fill borrowed slots with clones from the depleted cell
        c_idx = cell_sd[c]
        for (n, b) in enumerate(borrowed)
            src = c_idx[mod1(n, N_c)]
            ξ_src = droplets.ξ[src]
            # ξ_src == 2 && continue
            ξ_half = floor(Int, ξ_src / 2)
            droplets.ξ[src] = max(ξ_src - ξ_half, 1)
            droplets.ξ[b] = ξ_half
            droplets.X[b] = droplets.X[src]
            droplets.dry_r3[b] = droplets.dry_r3[src]
            droplets.cell_id[b] = c
            droplets.z_loc[b] = droplets.z_loc[src]
            droplets.w_prime[b] = FT(0)
        end

        # update snapshot so subsequent cells see the corrected donor count
        cell_sd[donor] = remaining
        append!(cell_sd[c], borrowed)
    end

    return nothing
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
        