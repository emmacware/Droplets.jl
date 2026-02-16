using RRTMGP



#these inputs have the same names as the scm.jl file so if you run that first you can
#call this function without it needing to be part of the timestepper

function stand_in_radiation_update_function!(grid,spatialsettings, diagnosticsettings, constants)
    #number of grid points and grid spacing
    Nz = grid.nz
    dz = grid.dz
    #Face (Nz+1) and center (Nz) heights
    z_faces = grid.faces_z
    z_centers = grid.centers_z

    #Profiles
    θ, qv, T, P, ρ = grid.states.θ, grid.states.qv, grid.states.T, grid.states.P, grid.states.ρ
    # plot_env_profiles(grid) will plot these

    #Liquid water diagnostics (cloud_, aerosol_, and rain_)
    cloud_LWC = grid.diagnostics.cloud_LWC
    cloud_effective_radius = grid.diagnostics.cloud_effective_radius

    return
end