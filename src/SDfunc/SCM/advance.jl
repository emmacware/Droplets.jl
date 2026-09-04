export single_column_timestep, update_theta!, update_density!, recycle_precipitation!, recycle_top_escape!, fill_grid_ranges!, keep_layer_filled!, radiation_function!

"""
radiation_function!
    Dispatches on the `radiation` dynamics flag to compute (or skip) radiative fluxes/heating
    for one timestep. The `DynOFF` no-op lives here since it needs no radiation-package
    dependency; the `DynON` method (RRTMGP-backed) is defined as `Droplets.radiation_function!`
    in the Examples radiation driver, so RRTMGP stays out of this package's dependencies.
"""
function radiation_function!(::DynOFF, grid::scm_eulerian_arrays{FT}, spatialsettings::spatial_settings{FT},
    diagnosticsettings::diagnostic_settings{FT}, constants::Constants{FT}, raddata, dt::FT, i::Int, n_rad::Int=1)::Nothing where FT<:AbstractFloat
    return nothing
end

"""
single_column_timestep

    Advance the SCM by one timestep, including condensation, coalescence, turbulence, advection, and diagnostics,
    if processes are turned on by scmsettings.

"""
function single_column_timestep(grid::scm_eulerian_arrays{FT}, dt::FT, droplets::droplet_attributes{FT},
    coagsettings::coag_settings{FT}, spatialsettings::spatial_settings{FT},
    condensationsettings::condensation_settings{FT}, coagdata::coagulation_run_spatial,
    conddata::condensation_data{FT}, raddata::rad,turbdata::turbulence_data{FT},
    diagnosticsettings::diagnostic_settings{FT},
    prescribed_w::Function, mpdatatmp::mpdata_tmp_1d, mpdatasettings::mpdata_settings_1d,
    constants::Constants{FT},
    scmsettings::scm_settings{FT},
    tkesettings::tke_settings{FT},
    absliq_r_interp,
    i::Int

) where {FT<:AbstractFloat, rad}
    nz = spatialsettings.Nz
    output_idx = Int(floor(i*dt/spatialsettings.dt_output) +1)
    n_output = spatialsettings.dt_output / dt

    #316.708 μs (0 allocations: 0 bytes)
    sd_fill_diagnostics(droplets, grid, spatialsettings, diagnosticsettings)
    # or i = 1
    if i % n_output == 0 || i == 1
        #11.666 μs (18 allocations: 10.67 KiB)
        scm_fill_diagnostic_output(grid,coagdata,conddata,raddata,spatialsettings,constants, output_idx, turbdata)
    end

    # 3.679 ms (2 allocations: 128 bytes) #REM off
    radiation_function!(scmsettings.radiation,grid, spatialsettings, diagnosticsettings, constants, raddata, dt, i, scmsettings.n_rad)

    T_from_theta!(grid.states.T_tmp,grid.states.θ, grid.states.P, grid.states.qv,constants)
    for _ in 1:scmsettings.n_cond
        #911.292 μs (48 allocations: 20.97 KiB)
        condensation_time_step_spatial!(scmsettings.condensation,droplets, grid.states, nz, dt/scmsettings.n_cond, conddata, constants, condensationsettings, spatialsettings, raddata,scmsettings,absliq_r_interp, i*dt)
    end




    coalescence_timestep!(scmsettings.coalescence,scmsettings.coag_threading, scmsettings.scheme, droplets, coagdata, coagsettings, scmsettings.n_coag)

    #2.990 μs (12 allocations: 5.31 KiB)
    T_from_theta!(grid.states.T_tmp,grid.states.θ, grid.states.P, grid.states.qv,constants)
    #10.000 μs (9 allocations: 4.28 KiB)
    compute_ql_at_cell!.(grid.states, 1:nz,constants)

    grid.states.qt_tmp .= grid.states.qv + grid.states.ql_tmp
    grid.states.θl_tmp .= θl.(grid.states.P, grid.states.T_tmp, grid.states.ql_tmp, grid.states.qv, constants)


    # Environmental advection
    #  160.250 μs (1799 allocations: 779.50 KiB)
    turb_timestep!(scmsettings.turbulence,grid, tkesettings, constants, dt, scmsettings, turbdata)

    #(32.916 μs (2 allocations: 1.06 KiB)
    compute_ql_all_cells!(grid.states, constants)
    if tkesettings.thermo_variable isa ThetalQtVar
        # turbulent diffusion above ran on θl_tmp/qt_tmp (see diffuse_fields! in tke.jl);
        # reconstruct θ/qv from the diffused pair + the freshly recomputed ql
        grid.states.qv .= grid.states.qt_tmp .- grid.states.ql_tmp
        grid.states.θ .= θ_from_θl.(grid.states.P, grid.states.θl_tmp, grid.states.ql_tmp, grid.states.qv, constants)
    end
    #88.291 ns (1 allocation: 64 bytes)
    update_theta!(scmsettings.density_feedback, grid, constants)
    #(3.010 μs (1 allocation: 64 bytes))
    update_density!(scmsettings.density_feedback, grid, constants)

    #33.167 μs (24 allocations: 13.30 KiB)
    mpdata_scm!(scmsettings.advection, scmsettings.thermo_feedback, grid, dt, mpdatatmp, mpdatasettings, constants)
    # 197.083 μs (16 allocations: 704 bytes)
    update_droplet_positions!(scmsettings.motion,scmsettings.advection,scmsettings.settling,droplets, prescribed_w, dt, spatialsettings, scmsettings, i)


    #(28.000 μs (49 allocations: 928 bytes))
    recycle_precipitation!(scmsettings.recycling,droplets, grid, spatialsettings, diagnosticsettings,coagdata, constants, output_idx)
    recycle_top_escape!(scmsettings.top_escape, droplets, spatialsettings, grid, constants)
    #18.833 μs (53 allocations: 1.03 KiB)
    keep_layer_filled!(scmsettings.keep_grid_filled, droplets, grid, spatialsettings, diagnosticsettings, constants, min_count = Int(floor(1.3 * coagsettings.Ns/spatialsettings.Nz)))
    #  12.916 μs (0 allocations: 0 bytes)
    fill_grid_ranges!(droplets)

return nothing
end



update_theta!(::DynON, grid, constants) = nothing
update_theta!(::DynOFF, grid, constants) = nothing

update_density!(::DynON, grid, constants) =
ρ_calc_θ!(grid.states.ρ,grid.states.P,grid.states.θ,grid.states.qv,constants)

# KiD: ρ_dry/P_dry were captured once at init and never touched;
# total ρ/P are re-diagnosed from the currently evolving
# qv (calc_ρ_from_ρ_dry / calc_P_from_P_dry)
function update_density!(::DynOFF, grid, constants)
states = grid.states
states.ρ .= calc_ρ_from_ρ_dry.(states.ρ_dry, states.qv)
states.P .= calc_P_from_P_dry.(states.P_dry, states.qv, Ref(constants))
return nothing
end

"""
recycle_precipitation
    Recycle precipitating droplets (X > diagnosticsettings.cloud_rain_cuttoff) that have fallen below the bottom of the grid
    back to the top of the grid, and add their mass to the surface precipitation diagnostic.
"""
function recycle_precipitation!(::DynOFF,droplets, grid, spatialsettings, diagnosticsettings, coagdata,
constants,output_step)
end

function recycle_precipitation!(::DynON, droplets, grid, spatialsettings, diagnosticsettings, coagdata, constants, output_step)
    FT = eltype(droplets.X)
    precip_mass = zero(FT)

    # grid_range may be stale after update_droplet_positions! moved droplets between
    # cells; refresh before using it to find drop_to_split
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

function recycle_top_escape!(::DynOFF, droplets, spatialsettings, grid, constants)
end
"""
recycle_top_escape
    Recycle droplets that have escaped the top of the grid back to the bottom of the grid.
    Setup for KiD
"""
function recycle_top_escape!(::DynON, droplets, spatialsettings, grid, constants)
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

"""
fill_grid_ranges!
    Fill the grid_range vector of UnitRanges for each grid cell, based on the current droplets.I and droplets.cell_id.
    This is used to efficiently access the droplets in each grid cell without having to filter the entire droplets.I array.
"""
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

function keep_layer_filled!(::DynOFF, droplets, grid, spatialsettings, diagnosticsettings, constants; min_count::Int = 100)
end
"""
keep_layer_filled!
    Ensure that each grid cell in the cloud layer of DYCOMS(between 400 m and the inversion height) has at least min_count superdroplets.
    If a cell has fewer than min_count superdroplets, it will borrow superdroplets from aerosol heavy cells with more than min_count superdroplets,
    superdroplets will be merged to free up slots, and the borrowed superdroplets will be cloned (and halved) from the depleted cell to maintain the 
    total number of superdroplets and real droplets.
    This is done to maintain a sufficient number of superdroplets for accurate representation of the droplet size distribution in each grid cell.
"""
function keep_layer_filled!(::DynON, droplets, grid, spatialsettings, diagnosticsettings, constants; min_count::Int = 100)
    FT = eltype(droplets.X)
    fill_grid_ranges!(droplets)  # ensure ranges are current
    nz = length(droplets.grid_range)

    dz = spatialsettings.z_grid_height
    z_inv_idx = findfirst(k -> grid.states.qv[k] <= FT(0.008), 1:length(droplets.grid_range)-1) - 1
    idx_400m = findfirst(k -> grid.centers_z[k] >= 400, 1:length(droplets.grid_range)-1)

    all(c -> length(droplets.grid_range[c]) >= min_count, idx_400m:z_inv_idx) && return nothing
    cell_sd = [collect(droplets.I[droplets.grid_range[g]]) for g in 1:nz]


    for c in z_inv_idx:-1:idx_400m
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
        for g in z_inv_idx+3:nz
            g == c && continue
            lg = length(cell_sd[g])
            if lg > donor_count
                donor_count = lg
                donor = g
            end
        end
        if donor == 0
            # println(N_needed)
            break
        end

        # println("Donor cell: ", donor, " with count: ", donor_count)

        # Check if donor is valid
        if donor == 0 || donor_count <= 40
            println(N_needed)
            # println("No valid donor found or donor count too low.")
            N_needed = 5
            # break
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
