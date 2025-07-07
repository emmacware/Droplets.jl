export mpdata_step!, flux_horizontal!, flux_vertical!, donor_cell_pass!, compute_antidiffusive_velocity!

limit(bc::Periodic,i,N) = i > N ? i + 1 - N : i < 1 ? i + N : i
limit(bc::NoFlux,i,N) = i > N ? N + 1 - (i - N) : i < 1 ? i + (1 + abs(i)) : i
flux(ψL, ψR, GC) = (max(GC, 0) * ψL + min(GC, 0) * ψR)

function flux_horizontal!(ϕ,ϕ_tmp,GCx;bc::BoundaryCondition=Periodic())
    nx, ny = size(ϕ)
    for i in 1:nx, j in 1:ny
        ip1,im1 = limit(bc,i+1, nx), limit(bc,i-1, nx)
        ϕ_tmp[i, j] -= flux(ϕ[i, j], ϕ[ip1, j], GCx[i+1, j])
        ϕ_tmp[i, j] += flux(ϕ[im1, j], ϕ[i, j], GCx[i, j])
    end
end

function flux_vertical!(ϕ,ϕ_tmp,GCy; bc::BoundaryCondition=Periodic())
    nx, ny = size(ϕ)
    for i in 1:nx, j in 1:ny   
        jp1,jm1 = limit(bc,j+1, ny), limit(bc,j-1, ny)
        ϕ_tmp[i, j] -= flux(ϕ[i, j], ϕ[i, jp1], GCy[i, j+1])
        ϕ_tmp[i, j] += flux(ϕ[i, jm1], ϕ[i, j], GCy[i, j])
    end
end

function donor_cell_pass!(ϕ,tmp; hbc::BoundaryCondition=Periodic(), vbc::BoundaryCondition=Periodic())
    tmp.ϕ .= 0

    flux_horizontal!(ϕ,tmp.ϕ,tmp.GCx_step,bc=hbc)
    flux_vertical!(ϕ,tmp.ϕ,tmp.GCy_step, bc=vbc)
    ϕ .+= tmp.ϕ
end


function u_corrective!(ϕ,tmp; f = 0.5,ϕ_eps=1e-20,hbc::BoundaryCondition=Periodic())
    GCy = tmp.GCy_step
    nx, ny = size(ϕ)
    for i in 1:nx, j in 2:ny-1

        ip, im = limit(hbc,i + 1, nx), limit(hbc,i - 1, nx)

        # A: avg
        A = (ϕ[ip,j] - ϕ[i,j]) / (ϕ[ip,j] + ϕ[i,j] + ϕ_eps)

        # B: (4-point stencil)
        jp1, jm1 = j+1, j-1
        B_num = (ϕ[ip,jp1] + ϕ[i,jp1]) - (ϕ[ip,jm1] + ϕ[i,jm1])
        B_den = 2 *(ϕ[ip,jp1] + ϕ[i,jp1] + ϕ[ip,jm1] + ϕ[i,jm1]) + ϕ_eps
        B = B_num / B_den

        U = tmp.GCx_step[ip,j]
        V = 0.25 * (GCy[i,j] + GCy[i,j+1] + GCy[ip,j] + GCy[ip,j+1]) 

        tmp.GCx_tmp[ip,j] = abs(U)*(1 - abs(U)) * A - 2*f*U*V*B
    end
end

function v_corrective!(ϕ,tmp; f = 0.5,ϕ_eps=1e-20, vbc::BoundaryCondition=Periodic())
    GCx = tmp.GCx_step
    nx, ny = size(ϕ)
    for i in 2:nx-1, j in 1:ny

        jp,jm = limit(vbc,j + 1, ny), limit(vbc,j - 1, ny)

        # B: avg
        B = (ϕ[i,jp] - ϕ[i,j]) / (ϕ[i,jp] + ϕ[i,j] + ϕ_eps)

        # A: (4-point stencil)
        ip1, im1 = i+1, i-1
        A_num = (ϕ[ip1,jp] + ϕ[ip1,j]) - (ϕ[im1,jp] + ϕ[im1,j])
        A_den = 2 * (ϕ[ip1,jp] + ϕ[ip1,j] + ϕ[im1,jp] + ϕ[im1,j]) + ϕ_eps
        A = A_num / A_den

        V = tmp.GCy_step[i,jp]
        U = 0.25 * (GCx[i,j] + GCx[ip1,j] + GCx[i,jp] + GCx[ip1,jp]) 

        tmp.GCy_tmp[i,jp] = abs(V)*(1 - abs(V)) * B - 2*(1 - f)*U*V*A
    end
end

function compute_antidiffusive_velocity!(ϕ, tmp; f = 0.5,hbc::BoundaryCondition = Periodic(), vbc::BoundaryCondition = Periodic())

    tmp.GCx_tmp .= 0
    tmp.GCy_tmp .= 0  

    u_corrective!(ϕ,tmp,hbc=hbc)
    v_corrective!(ϕ,tmp,vbc=vbc)

    tmp.GCx_step .= tmp.GCx_tmp
    tmp.GCy_step .= tmp.GCy_tmp
end

function mpdata_step!(ϕ_stage::Matrix{Float64}, GCx::Matrix{Float64}, GCy::Matrix{Float64}, tmp::mpdata_tmp,settings::mpdata_settings) #current 201.875 μs (0 allocations: 0 bytes)
    n_corr = settings.n_corr
    hbc::BoundaryCondition = settings.horizontal_boundary_condition
    vbc::BoundaryCondition = settings.vertical_boundary_condition

    tmp.GCx_step .= GCx.+0
    tmp.GCy_step .= GCy.+0

    donor_cell_pass!(ϕ_stage,tmp,hbc=hbc, vbc=vbc)

    for _ in 2:n_corr
        compute_antidiffusive_velocity!(ϕ_stage,tmp,hbc=hbc, vbc=vbc)
        donor_cell_pass!(ϕ_stage,tmp, hbc=hbc, vbc=vbc)
    end

end

function mpdata_step!(ϕ_stage::Tuple{Matrix{Float64}}, GCx::Matrix{Float64}, GCy::Matrix{Float64}, tmp::mpdata_tmp,settings::mpdata_settings) #current 201.875 μs (0 allocations: 0 bytes)
    n_corr = settings.n_corr
    hbc::BoundaryCondition = settings.horizontal_boundary_condition
    vbc::BoundaryCondition = settings.vertical_boundary_condition

    tmp.GCx_step .= GCx.+0
    tmp.GCy_step .= GCy.+0

    broadcast(ϕ -> donor_cell_pass!(ϕ, tmp, hbc=hbc, vbc=vbc), ϕ_stage)

    for _ in 2:n_corr
        broadcast(ϕ -> begin
            compute_antidiffusive_velocity!(ϕ, tmp, hbc=hbc, vbc=vbc)
            donor_cell_pass!(ϕ, tmp, hbc=hbc, vbc=vbc)
        end, ϕ_stage)
    end

end

# # Average to center
# function avg_u(u)
#     @views 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
# end
# function avg_v(v)
#     @views 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])
# end