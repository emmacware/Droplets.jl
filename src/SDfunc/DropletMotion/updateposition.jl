############################################################################################
# This file contains functions to update the position of the droplets in the domain 
############################################################################################

export update_position!, update_droplet_positions!

# ############################################################################################

# limit(a, N_x) = a > N_x ? a - N_x : a < 1 ? a + N_x : a
# limity(a, N_y) = a > N_y ? N_y : a < 1 ? 1 : a
# limitvface(a, N_y) = a > N_y +1 ? N_y +1 : a < 1 ? 1 : a


function update_droplet_positions!(droplets, w_function, Δt,spatialsettings)
    droplet_flow = w_function.(droplets.z_loc) .- v_term.(volume_to_radius.(droplets.X))
    droplets.z_loc .+= droplet_flow .* Δt
    droplets.cell_id .= Int.(floor.(droplets.z_loc ./ spatialsettings.z_grid_height) .+ 1)
    # need to handle droplets that have precipitated out of the domain
end

function update_gridflow_position!(droplets::droplet_attributes_1d,dt,grid_vels,spatial_settings)
    (;upper_face_vel, lower_face_vel) = grid_vels
    for droplet in droplets
        droplet.z_loc_in_cell += dt .* (
            (1 - droplet.z_loc_in_cell/ spatial_settings.z_grid_height ).* lower_face_vel .+
            (droplet.z_loc_in_cell / spatial_settings.z_grid_height) .* upper_face_vel
        )
    end
end

function update_gridflow_position!(droplets::droplet_attributes_2d,dt,grid_vels)
    (;upper_face_vel, lower_face_vel, right_face_vel, left_face_vel) = grid_vels
    for droplet in droplets
        droplet.z_loc_in_cell += dt .* (
            (1 - droplet.z_loc_in_cell / spatial_settings.z_grid_height) .* lower_face_vel .+
            (droplet.z_loc_in_cell / spatial_settings.z_grid_height) .* upper_face_vel
        )
        droplet.x_loc_in_cell += dt .* (
            (1 - droplet.x_loc_in_cell / spatial_settings.x_grid_width) .* left_face_vel .+
            (droplet.x_loc_in_cell / spatial_settings.x_grid_width) .* right_face_vel
        )
    end
end

function update_position!(droplets::Union{droplet_attributes_1d,droplet_attributes_2d},dt, spatial_settings,grid_vels)
    update_gridflow_position!(droplets,dt,grid_vels)

    if spatial_settings.settling == true
        droplets.z_loc_in_cell += dt .* terminal_velocity.(droplets.R)
    end
    #test for grid cell id changes.. or just recalculate?
    droplets.cell_id .= find_grid_cell_id.(droplets,spatial_settings)
end

function velocity_tendencies!(du,droplets::droplet_attributes_1d,spatial_settings)
    du.z_loc_in_cell .= find_grid_velocity_z.(droplets)
    if spatial_settings.settling == true
        du.z_loc_in_cell .+= terminal_velocity.(volume_to_radius(droplets.X))
    end
end
function velocity_tendencies!(du,droplets::droplet_attributes_2d,spatial_settings)
    du.z_loc_in_cell, du.x_loc_in_cell .= find_grid_velocity.(droplets)
    if spatial_settings.settling == true
        du.z_loc_in_cell .+= terminal_velocity.(volume_to_radius(droplets.X))
    end
end
    



