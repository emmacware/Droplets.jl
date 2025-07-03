############################################################################################
# This file contains functions to update the position of the droplets in the domain 
############################################################################################

export update_position!

# Spatial Settings Struct
Base.@kwdef struct spatial_settings{FT<:AbstractFloat}
    Nx::Int = 30 
    Nz::Int = 30 
    num_grids ::Int = Nx * Nz
    x_domain::FT = FT(1500.0) 
    z_domain::FT = FT(1500.0)
    z_grid_height::FT = z_domain / Nz 
    x_grid_width::FT = x_domain / Nx 
    periodic_boundaries_x::Bool = true 
    settling::Bool = true 
end
############################################################################################

limit(a, N_x) = a > N_x ? a - N_x : a < 1 ? a + N_x : a
limity(a, N_y) = a > N_y ? N_y : a < 1 ? 1 : a
limitvface(a, N_y) = a > N_y +1 ? N_y +1 : a < 1 ? 1 : a

function add_to_grid!(droplet, grid_index, grid_dict)
    if haskey(grid_dict, grid_index)
        push!(grid_dict[grid_index], droplet)
    else
        grid_dict[grid_index] = [droplet]
    end
end

# Function to move a superdroplet to a new grid
function move_to_grid!(droplet, old_grid_index, new_grid_index, grid_dict)
    # Remove droplet from old grid
    filter!(d -> d !== droplet, grid_dict[old_grid_index])

    # Add droplet to new grid
    add_to_grid!(droplet, new_grid_index, grid_dict)
end


function update_position!(droplets,Nx,Ny,dt,ρu,ρv,ρ, grid_dict,gridbox)
    moved_droplets = Set()

    for i in 1:Nx
        for j in 1:Ny
            if isempty(grid_dict[i,j])
                continue
            else
                n = limitvface(j+1,Ny)
                s = j
                e = limit(i+1,Nx)
                w = i
                for droplet in grid_dict[i,j]
                    # print(droplet.loc[1])
                    if droplet in moved_droplets
                        continue
                    else
                        droplet.loc[1] += dt*((droplet.loc[1]-i*Δx)*ρu[w,j]+ (1-(droplet.loc[1]-i*Δx))*ρu[e,j])/(2*ρ[i,j])
                        droplet.loc[2] += dt*(((droplet.loc[1]-i*Δx)*ρv[i,n]+ (1-(droplet.loc[2]-i*Δy))*ρv[i,s])/(2*ρ[i,j])-terminal_v(droplet.R))
                        push!(moved_droplets, droplet)
                        # print(droplet.loc[1])

                        if droplet.loc[1] <= gridbox[i,j][1] || droplet.loc[1] >= gridbox[i,j][2] 
                            if droplet.loc[1] < 0
                                droplet.loc[1] = droplet.loc[1]+Nx*Δx
                            elseif droplet.loc[1] > Nx*Δx
                                droplet.loc[1] = droplet.loc[1]-Nx*Δx
                            end
                            move_to_grid!(droplet, (i,j), (ceil(droplet.loc[1]/Δx),j), grid_dict)
                        end

                        if droplet.loc[2] <= gridbox[i,j][3]
                            if droplet.loc[2] <= 0
                                droplet.loc[2] = 0 #gridbox[i,j][3]
                                move_to_grid!(droplet, (i,j), (0,0), grid_dict)
                            else
                                move_to_grid!(droplet, (i,j), (i,(ceil(droplet.loc[2]/Δy))), grid_dict)
                            end
                        elseif droplet.loc[2]  >= gridbox[i,j][4] 
                            if i == Ny
                                droplet.loc[2] = gridbox[i,Ny][4]
                            else
                            move_to_grid!(droplet, (i,j), (i,(ceil(droplet.loc[2]/Δy))), grid_dict)
                            end
                        end
                    end
                end
            end
        end
    end
    return droplets, grid_dict
        
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

function find_grid_cell_id(droplets::droplet_attributes_1d,spatial_settings)
    for i in eachindex(droplets.z_loc_in_cell)
        if droplets.z_loc_in_cell[i] < 0
            droplets.cell_id[i] -= 1
            if droplets.cell_id != 0
                droplets.z_loc_in_cell[i] += spatial_settings.z_grid_height
            end
        elseif droplets.z_loc_in_cell[i] > spatial_settings.z_grid_height
            droplets.cell_id[i] += 1 
            if droplets.cell_id[i] != spatial_settings.Nz
                droplets.z_loc_in_cell[i] -= spatial_settings.z_grid_height
            end
        end
    end
end
    
function find_grid_cell_id(droplets::droplet_attributes_2d,spatial_settings)
    for i in eachindex(droplets.z_loc_in_cell)
        if droplets.z_loc_in_cell[i] < 0
            droplets.cell_id[i] -= Nz
            if droplets.cell_id <= 0
                droplets.cell_id = 0
            else
                droplets.z_loc_in_cell[i] += spatial_settings.z_grid_height
            end
        elseif droplets.z_loc_in_cell[i] > spatial_settings.z_grid_height
            droplets.cell_id[i] += Nz 
            if droplets.cell_id[i] >= spatial_settings.Nz * spatial_settings.Nx
                droplets.cell_id[i] = spatial_settings.Nz * spatial_settings.Ny + 1
            else
                droplets.z_loc_in_cell[i] -= spatial_settings.z_grid_height
            end
        end
    end

    for i in eachindex(droplets.x_loc_in_cell)
        if droplets.x_loc_in_cell[i] < 0
            if rem(droplets.cell_id[i], spatial_settings.Nx) == 1
                droplets.cell_id[i] += spatial_settings.Nx - 1
            else
                droplets.cell_id[i] -= 1
            end
        elseif droplets.x_loc_in_cell[i] > spatial_settings.x_grid_width
            if rem(droplets.cell_id[i], spatial_settings.Nx) == 0
                droplets.cell_id[i] -= (spatial_settings.Nx - 1)
            else
                droplets.cell_id[i] += 1
            end
        end
    end
end