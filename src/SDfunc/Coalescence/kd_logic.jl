#---------------------------------------------------------
# KinematicDriver.jl compatibility logic
#---------------------------------------------------------
using StaticArrays

export coalescence_timestep!, static_droplet_attributes
export KiD, spatial, coalescence_timestep!, size_thresh_separate_droplets

struct KiD end #<: scheme_type end
struct spatial end #<: scheme_type



function coalescence_timestep!(run::Serial,scheme::KiD, ξFT::SVector,X::SVector,
    coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    
    Ns::Int = 0  # Declare once with a default value
    first_unhealthy_idx = findfirst(iszero, ξFT)
    if first_unhealthy_idx === 1
        return (; SD_Vol = X, SD_Mult = ξFT)
    elseif first_unhealthy_idx === nothing
        Ns = length(ξFT)
    else
        Ns = first_unhealthy_idx - 1
    end
    
    droplets = droplet_attributes{FT}(Int.(floor.(Vector(ξFT))), Vector(X))
    
    I = shuffle(1:Ns)
    L = [(I[l-1], I[l]) for l in 2:2:Ns]
    scale::FT = Ns * (Ns - 1) / 2 / (Ns / 2)

    compute_pαdt!(L, droplets,coag_data,settings.kernel,settings)

    rand!(coag_data.ϕ)

    test_pairs!(run,L,droplets,coag_data)
    ξFT = SVector{length(ξFT), FT}(FT.(droplets.ξ))
    X = SVector{length(ξFT), FT}(droplets.X)
    return (;SD_Vol = X,SD_Mult = ξFT)
end 


function size_thresh_separate_droplets(SD_Vol,SD_Mult,ρd)
    FT = eltype(SD_Vol)
    # Define the threshold for separating droplets
    threshold = 40*1e-6 # meters
    threshold_volume = (4/3) * π * (threshold^3) # m^3

    # Separate the droplets based on the threshold
    N_liq = sum(SD_Mult[SD_Vol .< threshold_volume])
    N_rai = sum(SD_Mult[SD_Vol .>= threshold_volume])

    q_liq = 1000*sum(SD_Vol[SD_Vol .< threshold_volume].* SD_Mult[SD_Vol .< threshold_volume])/ρd
    q_rai = 1000*sum(SD_Vol[SD_Vol .>= threshold_volume].* SD_Mult[SD_Vol .>= threshold_volume])/ρd

    return (N_liq=FT(N_liq), N_rai=FT(N_rai), q_liq=FT(q_liq), q_rai=FT(q_rai))
end




function coalescence_timestep!(run::Serial,scheme::spatial,spatial_settings, droplets)
    
    filled = 0 
    L = [(0, 0) for l in 2:2:Ns]
    shuffle!(coag_data.I)
    for cell in 1:spatial_settings.Nx*spatial_settings.Ny
        pairs_in_cell = findall(x ->droplets.cell_id[x] == cell, I) #replace with dict?
        working_pairs = div(length(pairs_in_cell), 2)
        L[filled:filled+working_pairs] = [(pairs_in_cell[l-1], pairs_in_cell[l]) for l in 2:2:length(pairs_in_cell)]
        coag_data.scale[filled:filled+working_pairs] = working_pairs * (working_pairs - 1) / 2 / (Ns / 2)
        filled += working_pairs
    end
    L = L[1:filled]

    compute_pαdt!(L, droplets,coag_data,settings.kernel,settings)
    rand!(coag_data.ϕ)
    test_pairs!(run,L,droplets,coag_data)
end 

#---------------------------------------------------------

# droplets operator splittings 

function rhs_droplets!(du,Y,p,t)
    FT = eltype(Y)
    (;run_settings,spatial_settings,constants,droplets) = p
    (;R_air, cp_air, L_vap, ρ_air) = p
    velocity_field = p
    X = Y.X

    for i in 1:spatial_settings.num_grids
        Senv = Y_CompVec.Senv[i]
        T = Y_CompVec.T[i]
        qv = Y_CompVec.qv[i]

        in_cell = findall(x -> droplets.cell_id[x] == i, 1:length(X))
        if isempty(in_cell)
            continue
        end
        # Condensation/evaporation
        dX = dX_droplets!(X[in_cell],droplets.dry_r3[in_cell], run_settings.kappa, qv, T, velocity_field.P, t)
        dql = sum(dX.*droplets.ξ[in_cell].* constants.ρl / ρ_air) 
        dqv = -dql

        dSenv = - (1/qv + (L_vap^2 / constants.Rv / T^2 / cp_air)) * Senv * dql

        dT = L_vap / cp_air * dql

        du.X[in_cell] .= dX
        du.Y_CompVec.Senv[i] .= dSenv
        du.Y_CompVec.T[i] .= dT
        du.Y_CompVec.qv[i] .= dqv

    end
end

# callback function no solver? maybe leave condensation separate, just motion + coalescence
function droplets_operater_splitting!(droplets, p, dt)
    (;run_settings, spatial_settings) = p

    update_position!(droplets::Union{droplet_attributes_1d,droplet_attributes_2d},dt, spatial_settings)

    coalescence_timestep!(run_settings.coag_threading, run_settings.scheme, droplets, run_settings.coagulation_data, run_settings.coagulation_settings)

end


function combined_timestepper!()


    #step clima ode

    #take info from clima Y: air densities, velocity fields, 

    #step condensation ode

    #move droplets and coalesce

    #return all the data to og integrator

end