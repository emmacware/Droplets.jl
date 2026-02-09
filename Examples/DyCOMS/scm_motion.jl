
function update_droplet_positions(droplets, w_function, Δt)
    droplet_flow = w_function.(droplets.z_loc) .- v_term.(volume_to_radius.(droplets.X))
    droplets.z_loc .+= droplet_flow .* Δt
    droplets.cell_id .= floor.(droplets.z_loc ./ grid.dz) .+ 1
    
    return
end


function mpdata_scm(grid, Δt, tmp, run_settings,constants)

    GCz = grid.wind.w * Δt ./ grid.dz

    mpdata_step!(grid.states.qv, GCz,tmp,run_settings)
    mpdata_step!(grid.states.θ, GCz,tmp,run_settings)
    
    grid.states.T .= in_situ_temperature(grid.states.θ,grid.states.qv,grid.states.P,constants)
    grid.states.ρ .= grid.states.P ./ (constants.Rd .* grid.states.T .* (1 .+ 0.61 .* grid.states.qv))
    return
end