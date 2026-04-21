export mpdata_step!, flux_horizontal!, flux_vertical!, donor_cell_pass!, compute_antidiffusive_velocity!
export antiosc!, flux_faces, beta_function, beta_function_1d, find_extrema!, mpdata_scm!

limit(bc::Periodic,i,N) = i > N ? i - N : i < 1 ? i + N : i
limit(bc::NoFlux,i,N) = i > N ? N + 1 - (i - N) : i < 1 ? i + (1 + abs(i)) : i
flux(ψL, ψR, GC;infinite_gauge::Bool=false) = infinite_gauge ? GC : (max(GC, 0) * ψL + min(GC, 0) * ψR)

function flux_horizontal!(ϕ,ϕ_tmp,GCx;bc::BoundaryCondition=Periodic(), infinite_gauge::Bool=false)
    nx, ny = size(ϕ)
    for i in 1:nx, j in 1:ny
        ip1,im1 = limit(bc,i+1, nx), limit(bc,i-1, nx)
        ϕ_tmp[i, j] -= flux(ϕ[i, j], ϕ[ip1, j], GCx[i+1, j],infinite_gauge=infinite_gauge)
        ϕ_tmp[i, j] += flux(ϕ[im1, j], ϕ[i, j], GCx[i, j],infinite_gauge=infinite_gauge)
    end
end

function flux_vertical!(ϕ,ϕ_tmp,GCy; bc::BoundaryCondition=Periodic(), infinite_gauge::Bool=false)
    nx, ny = size(ϕ)
    for i in 1:nx, j in 1:ny   
        jp1,jm1 = limit(bc,j+1, ny), limit(bc,j-1, ny)
        ϕ_tmp[i, j] -= flux(ϕ[i, j], ϕ[i, jp1], GCy[i, j+1],infinite_gauge=infinite_gauge)
        ϕ_tmp[i, j] += flux(ϕ[i, jm1], ϕ[i, j], GCy[i, j],infinite_gauge=infinite_gauge)
    end

end


function donor_cell_pass!(ϕ,tmp::mpdata_tmp_1d; vbc::BoundaryCondition=Periodic(), infinite_gauge::Bool=false,topcellreservoir=false)
    tmp.ϕ .= 0

    nz = length(ϕ)
    for k in 1:nz
        kp,km = limit(vbc,k+1, nz), limit(vbc,k-1, nz)
        tmp.ϕ[k] -= flux(ϕ[k], ϕ[kp], tmp.GCz_step[k+1],infinite_gauge=infinite_gauge)
        tmp.ϕ[k] += flux(ϕ[km], ϕ[k], tmp.GCz_step[k],infinite_gauge=infinite_gauge)
    end
    if topcellreservoir
        tmp.ϕ[end] = 0.0
        # tmp.ϕ[end] += flux(ϕ[end-1], ϕ[end], tmp.GCz_step[end-1],infinite_gauge=infinite_gauge)
    end

    ϕ .+= tmp.ϕ
end

function donor_cell_pass!(ϕ,tmp::mpdata_tmp; hbc::BoundaryCondition=Periodic(), vbc::BoundaryCondition=Periodic(), infinite_gauge::Bool=false)
    tmp.ϕ .= 0
    
    flux_horizontal!(ϕ,tmp.ϕ,tmp.GCx_step, bc=vbc,infinite_gauge=infinite_gauge)
    flux_vertical!(ϕ,tmp.ϕ,tmp.GCy_step, bc=vbc,infinite_gauge=infinite_gauge)

    ϕ .+= tmp.ϕ
end


function u_corrective!(ϕ,tmp,settings; f = 0.5,ϕ_eps=1e-20,hbc::BoundaryCondition=Periodic(),vbc::BoundaryCondition=Periodic())
    GCx_inside = tmp.GCx_step .+ 0
    GCy_inside = tmp.GCy_step .+ 0
    nx, ny = size(ϕ)

    for i in 1:nx, j in 1:ny

        ip = limit(hbc,i + 1, nx)
        jp1, jm1 = limit(vbc,j + 1, ny), limit(vbc,j - 1, ny)

        # A: avg
        A_num =(ϕ[ip,j] - ϕ[i,j])
        A_den = settings.infinite_gauge ? 2 : (ϕ[ip,j] + ϕ[i,j] + ϕ_eps)
        A = A_num / A_den

        # B: (4-point stencil)
        B_num = (ϕ[ip,jp1] + ϕ[i,jp1]) - (ϕ[ip,jm1] + ϕ[i,jm1])
        B_den = settings.infinite_gauge ? 8 : (2 * (ϕ[ip,jp1] + ϕ[i,jp1] + ϕ[ip,jm1] + ϕ[i,jm1]) + ϕ_eps )
        B = B_num /B_den

        U = GCx_inside[i+1,j]
        V = 0.25 * (GCy_inside[i,j] + GCy_inside[i,j+1] + GCy_inside[ip,j] + GCy_inside[ip,j+1]) 

        tmp.GCx_tmp[i+1,j] = (abs(U)*(1 - abs(U)) * A - 2*f*U*V*B)

    end

    if hbc isa Periodic
        tmp.GCx_tmp[1,:] .= tmp.GCx_tmp[end,:]
    end

end

function v_corrective!(ϕ,tmp,settings; f = 0.5,ϕ_eps=1e-20, vbc::BoundaryCondition=Periodic(),hbc::BoundaryCondition=Periodic())
    GCx_inside = tmp.GCx_step .+ 0
    GCy_inside = tmp.GCy_step .+ 0
    nx, ny = size(ϕ)

    for i in 1:nx, j in 1:ny

        jp= limit(vbc,j + 1, ny)
        ip1, im1 = limit(hbc,i + 1, nx), limit(hbc,i - 1, nx)

        # B: avg
        B_num = (ϕ[i,jp] - ϕ[i,j])
        B_den = settings.infinite_gauge ? 2 : (ϕ[i,jp] + ϕ[i,j] + ϕ_eps)
        B = B_num / B_den

        # A: (4-point stencil)
        A_num = (ϕ[ip1,jp] + ϕ[ip1,j]) - (ϕ[im1,jp] + ϕ[im1,j])
        A_den = settings.infinite_gauge ? 8 : (2 * (ϕ[ip1,jp] + ϕ[ip1,j] + ϕ[im1,jp] + ϕ[im1,j]) + ϕ_eps)
        A = A_num / A_den

        V = GCy_inside[i,j+1]
        U = 0.25 * (GCx_inside[i,j] + GCx_inside[i+1,j] + GCx_inside[i,jp] + GCx_inside[i+1,jp]) 

        tmp.GCy_tmp[i,j+1] = (abs(V)*(1 - abs(V)) * B - 2*(1 - f)*U*V*A)

    end

    if vbc isa Periodic
        tmp.GCy_tmp[:,1] .= tmp.GCy_tmp[:,end]
    end
end

function corrective_vel_1d(ϕ,tmp,settings; f = 0.5,ϕ_eps=1e-20, vbc::BoundaryCondition=Periodic())
    GCz_inside = tmp.GCz_step .+ 0
    nz = length(ϕ)

    for k_ in 0:nz-1

        kp = limit(vbc,k_ + 1, nz)
        k = limit(vbc,k_, nz)

        B_num =(ϕ[kp] - ϕ[k])
        B_den = settings.infinite_gauge ? 2 : (ϕ[kp] + ϕ[k] + ϕ_eps)
        B = B_num / B_den

        tmp.GCz_tmp[kp] = abs(GCz_inside[kp])*(1 - abs(GCz_inside[kp])) * B 
    end
    if vbc isa Periodic
        tmp.GCz_tmp[1] = tmp.GCz_tmp[end]
    end
end

function beta_function(ϕ,tmp,i_idx,j_idx,GCxph,GCxmh,GCyph,GCymh,infinite_gauge;ϕ_eps=1e-20)
    im1,i,ip1 = i_idx
    jm1,j,jp1 = j_idx
    ψ = ϕ[i,j]
    ψL = ϕ[im1,j]
    ψR = ϕ[ip1,j]
    ψU = ϕ[i,jp1]
    ψD = ϕ[i,jm1]
    max_phi, min_phi = tmp.minmax.localmax, tmp.minmax.localmin

    β_in = (max(max_phi[i,j],ψR, ψ,ψL,ψU,ψD) - ψ)/ (
        max(flux(ψL, ψ, GCxmh,infinite_gauge=infinite_gauge),0) +
        -min(flux(ψ, ψR, GCxph,infinite_gauge=infinite_gauge),0) +
        max(flux(ψD, ψ, GCymh,infinite_gauge=infinite_gauge),0) +
        -min(flux(ψ, ψU, GCyph,infinite_gauge=infinite_gauge),0) + ϕ_eps)

    β_out = (ψ - min(min_phi[i,j],ψR, ψ,ψL,ψU,ψD ))/ (
        -min(flux(ψL, ψ, GCxmh,infinite_gauge=infinite_gauge),0) +
        max(flux(ψ, ψR, GCxph,infinite_gauge=infinite_gauge),0) +
        -min(flux(ψD, ψ, GCymh,infinite_gauge=infinite_gauge),0) +
        max(flux(ψ, ψU, GCyph,infinite_gauge=infinite_gauge),0) + ϕ_eps)

    #limit beta to 0-1
    β_in = max(0.0, min(1.0, β_in))
    β_out = max(0.0, min(1.0, β_out))
    return β_in, β_out
end





# function flux_faces(ϕ,GCx,GCy;hbc::BoundaryCondition=Periodic(),vbc::BoundaryCondition=Periodic(), infinite_gauge::Bool=false)
#     nx, ny = size(ϕ)
#     flux_h,flux_v, = zeros(nx+1, ny), zeros(nx, ny+1)
#     for i in 1:nx, j in 1:ny
#         ip1 = limit(hbc,i+1, nx)
#         jp1 = limit(vbc,j+1, ny)

#         flux_h[i+1, j] = flux(ϕ[i, j], ϕ[ip1, j], GCx[i+1, j],infinite_gauge=infinite_gauge)
#         flux_v[i, j+1] = flux(ϕ[i, j], ϕ[i, jp1], GCy[i, j+1],infinite_gauge=infinite_gauge)
#     end
#     if hbc isa Periodic
#         flux_h[1,:] .= flux_h[end,:]
#     end
#     if vbc isa Periodic
#         flux_v[:,1] .= flux_v[:,end]
#     end
#     return flux_h, flux_v
# end



function antiosc!(ϕ,tmp,settings; hbc::BoundaryCondition=Periodic(), vbc::BoundaryCondition=Periodic())
    nx, ny = size(ϕ)
    
    tmp.GCx_tmp .= 0
    tmp.GCy_tmp .= 0 
    
    # # Compute fluxes
    # flux_h, flux_v = flux_faces(ϕ, antidiffusive_GCx, antidiffusive_GCy; 
    #                            hbc=hbc, vbc=vbc, infinite_gauge=settings.infinite_gauge)

    for i in 1:nx, j in 1:ny 
            jp1,jm1,jp2 = limit(vbc,j + 1, ny), limit(vbc,j - 1, ny), limit(vbc,j + 2, ny)
            ip1, im1,ip2 = limit(hbc,i + 1, nx), limit(hbc,i - 1, nx), limit(hbc,i + 2, nx)
            #faces
            v_jp2 = limit(vbc,j + 2, ny)
            u_ip2 = limit(hbc,i + 2, nx)

        beta_in,beta_out = beta_function(ϕ,tmp, (im1,i,ip1), (jm1,j,jp1),
        tmp.GCx_step[i+1,j],tmp.GCx_step[i,j],
        tmp.GCy_step[i,j+1],tmp.GCy_step[i,j],
                        settings.infinite_gauge)

        betajp1_in,betajp1_out = beta_function(ϕ,tmp, (im1,i,ip1), (j,jp1,jp2),
        tmp.GCx_step[i+1,jp1],tmp.GCx_step[i,jp1],
        tmp.GCy_step[i,v_jp2],tmp.GCy_step[i,j+1],
                        settings.infinite_gauge)
        
        betaip1_in,betaip1_out = beta_function(ϕ,tmp,(i,ip1,ip2), (jm1,j,jp1),
                        tmp.GCx_step[u_ip2,j],tmp.GCx_step[i+1,j],
                        tmp.GCy_step[ip1,j+1],tmp.GCy_step[ip1,j],
                        settings.infinite_gauge)

        
        # Calculate scaling factors instead of new fluxes
        # if abs(tmp.GCy_step[i,j+1]) > 1e-20
            if tmp.GCy_step[i,j+1] >= 0
                scale_v = min(1, beta_out, betajp1_in)
            else
                scale_v = min(1, beta_in, betajp1_out)
            end
            tmp.GCy_tmp[i,j+1] = tmp.GCy_step[i,j+1] * scale_v
        # else
        #     tmp.GCy_tmp[i,j+1] = 0.0
        # end
        
        # if abs(tmp.GCx_step[i+1,j]) > 1e-20
            if tmp.GCx_step[i+1,j] >= 0
                scale_h = min(1, beta_out, betaip1_in)
            else
                scale_h = min(1, beta_in, betaip1_out)
            end
            tmp.GCx_tmp[i+1,j] = tmp.GCx_step[i+1,j] * scale_h
        # else
        #     tmp.GCx_tmp[i+1,j] = 0.0
        # end
    end

    # Apply periodic BCs
    if hbc isa Periodic
        tmp.GCx_tmp[1,:] .= tmp.GCx_tmp[end,:]
    end
    if vbc isa Periodic
        tmp.GCy_tmp[:,1] .= tmp.GCy_tmp[:,end]
    end
    tmp.GCx_step .= tmp.GCx_tmp
    tmp.GCy_step .= tmp.GCy_tmp
end

function compute_antidiffusive_velocity!(ϕ, tmp::mpdata_tmp,settings; f = 0.5,hbc::BoundaryCondition = Periodic(), vbc::BoundaryCondition = Periodic())

    tmp.GCx_tmp .= 0
    tmp.GCy_tmp .= 0  

    u_corrective!(ϕ,tmp,settings,f=f,hbc=hbc)
    v_corrective!(ϕ,tmp,settings,f=f,vbc=vbc)

    tmp.GCx_step .= tmp.GCx_tmp
    tmp.GCy_step .= tmp.GCy_tmp
end

function compute_antidiffusive_velocity!(ϕ, tmp::mpdata_tmp_1d,settings; f = 0.5,vbc::BoundaryCondition = Periodic())

    tmp.GCz_tmp .= 0

    corrective_vel_1d(ϕ,tmp,settings,f=f,vbc=vbc)

    tmp.GCz_step .= tmp.GCz_tmp
end

function find_extrema!(ϕ,tmp;hbc::BoundaryCondition=Periodic(), vbc::BoundaryCondition=Periodic())
    nx, ny = size(ϕ)
    for i in 1:nx, j in 1:ny
        im1, ip1 = limit(hbc,i - 1, nx), limit(hbc,i + 1, nx)
        jm1, jp1 = limit(vbc,j - 1, ny), limit(vbc,j + 1, ny)
        tmp.minmax.localmax[i,j] = maximum((ϕ[im1,j], ϕ[i,j], ϕ[ip1,j], ϕ[i,jp1], ϕ[i,jm1]))
        tmp.minmax.localmin[i,j] = minimum((ϕ[im1,j], ϕ[i,j], ϕ[ip1,j], ϕ[i,jp1], ϕ[i,jm1]))
    end
end

function find_extrema!(ϕ::Vector, tmp::mpdata_tmp_1d; vbc::BoundaryCondition=Periodic())
    nz = length(ϕ)
    for k in 1:nz
        km1, kp1 = limit(vbc, k-1, nz), limit(vbc, k+1, nz)
        tmp.minmax.localmax[k] = maximum((ϕ[km1], ϕ[k], ϕ[kp1]))
        tmp.minmax.localmin[k] = minimum((ϕ[km1], ϕ[k], ϕ[kp1]))
    end
end

function beta_function_1d(ϕ, tmp::mpdata_tmp_1d, k, km1, kp1, GCzmh, GCzph, infinite_gauge; ϕ_eps=1e-20)
    ψ = ϕ[k]; ψL = ϕ[km1]; ψR = ϕ[kp1]
    max_phi = tmp.minmax.localmax[k]
    min_phi = tmp.minmax.localmin[k]

    β_in = (max_phi - ψ) / (
         max(flux(ψL, ψ, GCzmh, infinite_gauge=infinite_gauge), 0) +
        -min(flux(ψ, ψR, GCzph, infinite_gauge=infinite_gauge), 0) + ϕ_eps)

    β_out = (ψ - min_phi) / (
        -min(flux(ψL, ψ, GCzmh, infinite_gauge=infinite_gauge), 0) +
         max(flux(ψ, ψR, GCzph, infinite_gauge=infinite_gauge), 0) + ϕ_eps)

    β_in  = max(0.0, min(1.0, β_in))
    β_out = max(0.0, min(1.0, β_out))
    return β_in, β_out
end

function antiosc!(ϕ::Vector, tmp::mpdata_tmp_1d, settings; vbc::BoundaryCondition=Periodic())
    nz = length(ϕ)
    tmp.GCz_tmp .= 0

    for k in 1:nz-1
        km1 = limit(vbc, k-1, nz)
        kp1 = k + 1
        kp2 = limit(vbc, k+2, nz)

        beta_k_in,   beta_k_out   = beta_function_1d(ϕ, tmp, k,   km1, kp1,
            tmp.GCz_step[k],   tmp.GCz_step[k+1], settings.infinite_gauge)
        beta_kp1_in, beta_kp1_out = beta_function_1d(ϕ, tmp, kp1, k,   kp2,
            tmp.GCz_step[k+1], tmp.GCz_step[k+2], settings.infinite_gauge)

        if tmp.GCz_step[k+1] >= 0
            scale = min(1.0, beta_k_out, beta_kp1_in)
        else
            scale = min(1.0, beta_k_in,  beta_kp1_out)
        end
        tmp.GCz_tmp[k+1] = tmp.GCz_step[k+1] * scale
    end

    if vbc isa Periodic
        tmp.GCz_tmp[1] = tmp.GCz_tmp[end]
    end

    tmp.GCz_step .= tmp.GCz_tmp
end

function mpdata_step!(ϕ_stage::Matrix{Float64}, GCx::Matrix{Float64}, GCy::Matrix{Float64}, tmp::mpdata_tmp,settings::mpdata_settings;f=0.5) #current 201.875 μs (0 allocations: 0 bytes)
    n_corr = settings.n_corr
    hbc::BoundaryCondition = settings.horizontal_boundary_condition
    vbc::BoundaryCondition = settings.vertical_boundary_condition

    tmp.GCx_step .= GCx.+0
    tmp.GCy_step .= GCy.+0
    find_extrema!(ϕ_stage,tmp)

    donor_cell_pass!(ϕ_stage,tmp,hbc=hbc, vbc=vbc, infinite_gauge=false)


    for _ in 2:n_corr
        compute_antidiffusive_velocity!(ϕ_stage,tmp,settings,f=f,hbc=hbc, vbc=vbc)
        if settings.nonoscillatory
            antiosc!(ϕ_stage,tmp,settings,hbc=hbc,vbc=vbc)
            find_extrema!(ϕ_stage,tmp)
        end

        donor_cell_pass!(ϕ_stage,tmp, hbc=hbc, vbc=vbc,infinite_gauge=settings.infinite_gauge)

    end

end

# function mpdata_step!(ϕ_stage::Tuple{Matrix{Float64}}, GCx::Matrix{Float64}, GCy::Matrix{Float64}, tmp::mpdata_tmp,settings::mpdata_settings) #current 201.875 μs (0 allocations: 0 bytes)
#     n_corr = settings.n_corr
#     hbc::BoundaryCondition = settings.horizontal_boundary_condition
#     vbc::BoundaryCondition = settings.vertical_boundary_condition

#     tmp.GCx_step .= GCx.+0
#     tmp.GCy_step .= GCy.+0

#     broadcast(ϕ -> donor_cell_pass!(ϕ, tmp, hbc=hbc, vbc=vbc), ϕ_stage)

#     for _ in 2:n_corr
#         broadcast(ϕ -> begin
#             compute_antidiffusive_velocity!(ϕ, tmp, hbc=hbc, vbc=vbc)

#             donor_cell_pass!(ϕ, tmp, hbc=hbc, vbc=vbc)
#         end, ϕ_stage)
#     end

# end

# # Average to center
# function avg_u(u)
#     @views 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
# end
# function avg_v(v)
#     @views 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])
# end

function mpdata_step!(ϕ_stage::Vector{Float64}, GCz::Vector{Float64}, tmp::mpdata_tmp_1d,settings::mpdata_settings_1d;f=0.5) 
    n_corr = settings.n_corr
    vbc::BoundaryCondition = settings.vertical_boundary_condition

    tmp.GCz_step .= GCz.+0
    find_extrema!(ϕ_stage,tmp,vbc=vbc)

    tcr = vbc isa NoFlux
    donor_cell_pass!(ϕ_stage,tmp,vbc=vbc, infinite_gauge=false, topcellreservoir=tcr)


    for _ in 2:n_corr
        compute_antidiffusive_velocity!(ϕ_stage,tmp,settings,f=f, vbc=vbc)
        if settings.nonoscillatory
            antiosc!(ϕ_stage,tmp,settings,vbc=vbc)
            find_extrema!(ϕ_stage,tmp,vbc=vbc)
        end

        donor_cell_pass!(ϕ_stage,tmp, vbc=vbc,infinite_gauge=settings.infinite_gauge, topcellreservoir=tcr)

    end
end

function mpdata_scm!(grid::scm_eulerian_arrays, Δt::FT, tmp::mpdata_tmp_1d, settings::mpdata_settings_1d,constants::Constants) where {FT<:AbstractFloat}

    GCz = grid.wind.w * Δt ./ grid.dz

    ρqv = grid.states.ρ .* grid.states.qv
    ρθ = grid.states.ρ .* grid.states.θ
    ρe = grid.states.ρ .* grid.states.e

    # mpdata_step!(grid.states.qv, GCz,tmp,settings)
    # mpdata_step!(grid.states.θ, GCz,tmp,settings)
    mpdata_step!(ρqv, GCz,tmp,settings)
    mpdata_step!(ρθ, GCz,tmp,settings)
    mpdata_step!(ρe, GCz,tmp,settings)
    mpdata_step!(grid.states.ρ, GCz,tmp,settings)
    # mpdata_step!(grid.states.e, GCz,tmp,settings)
    grid.states.qv .= ρqv ./ grid.states.ρ
    grid.states.θ .= ρθ ./ grid.states.ρ
    grid.states.e .= ρe ./ grid.states.ρ
    
    # grid.states.T .= T_from_theta(grid.states.θ,grid.states.P,constants)
    grid.states.ρ .= ρ_calc_θ(grid.states.P,grid.states.θ,grid.states.qv,constants)
    return
end
