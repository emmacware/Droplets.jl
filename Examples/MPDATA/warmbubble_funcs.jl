

matrix = zeros(10,10)
for i in 1:10, j in 1:10
    matrix[i,j] = i+j
end

# matrix_arrays = (lap_tmp = zeros((2,10,10)),
#                 tmp_uvw = zeros((2,10,10)))

# laplacian(matrix, false, matrix_arrays, 20.0, 20.0,halo=2)
# grad(matrix, 1,20.0,halo=2)

function grad(phi, axis,dl; halo::Int=1)
    size_dim1, size_dim2 = size(phi)
    xchng_pres_old!(phi, halo)
    if axis == 1
        gradright = ((halo + 1+1):(size_dim1 - halo+1),
            (halo + 1):(size_dim2 - halo))
        gradleft = ((halo):(size_dim1 - halo-1),
            (halo + 1):(size_dim2 - halo))
    elseif axis == 2
        gradright = ((halo + 1):(size_dim1 - halo), 
        (halo + 1+1):(size_dim2 - halo+1))
        gradleft = (
            (halo + 1):(size_dim1 - halo), 
            (halo):(size_dim2 - halo - 1)
        )
    else
        error("Axis must be 1 or 2")
    end
    return (phi[gradright...] .- phi[gradleft...])./ 2 ./ dl
end


function laplacian(Phi::AbstractMatrix, err_init::Bool, ARRAYS, dx::Real, dy::Real; halo::Int=1)
    Nx_tot, Ny_tot = size(Phi)
    i0, i1 = 1+halo, Nx_tot-halo
    j0, j1 = 1+halo, Ny_tot-halo
    interior_indices = CartesianIndices((i0:i1, j0:j1))

    # ensure Phi halos are up-to-date
    xchng_pres_old!(Phi, halo)


    calc_grad_phi(ARRAYS.lap_tmp, Phi, interior_indices, dx, dy,halo=halo)
    xchng_pres_old!(@view ARRAYS.lap_tmp[1,:,:])
    xchng_pres_old!(@view ARRAYS.lap_tmp[2,:,:])
    xchng_pres_old!(@view ARRAYS.tmp_uvw[1,:,:])
    xchng_pres_old!(@view ARRAYS.tmp_uvw[2,:,:])

    # optional subtraction
    if err_init
        ARRAYS.lap_tmp[1,:,:] -= ARRAYS.tmp_uvw[1,:,:]
        ARRAYS.lap_tmp[2,:,:] -= ARRAYS.tmp_uvw[2,:,:]
    end

    # update halos of lap_tmp
    xchng_pres_old!(@view(ARRAYS.lap_tmp[1,:,:]), halo)
    xchng_pres_old!(@view(ARRAYS.lap_tmp[2,:,:]), halo)

    lap = grad(ARRAYS.lap_tmp[1,:,:], 1,dx,halo=halo) + grad(ARRAYS.lap_tmp[2,:,:], 2,dy,halo=halo)

    return lap
end

function calc_grad_phi(arg, phi,idx_interior,dx,dy; halo::Int=1)
    arg[1,idx_interior] .= grad(phi, 1,dx,halo=halo)
    arg[2,idx_interior] .= grad(phi, 2,dy,halo=halo)
end

function update_rhs(ARRAYS,θ,g,θ_ref)
    ARRAYS.rhs_w .= g .* (θ .- θ_ref) ./ θ_ref
end


function pressure_solver_loop_body!(ARRAYS, converged,k_iters=4,error_tol=1e-7/1.5)
        # Setup temporary arrays
        tmp_den = ones(k_iters)
        alpha   = ones(k_iters)
        nx, ny = size(ARRAYS.lap_err)
        #create the slices for interior points
        idx_interior = CartesianIndices((2:nx+1, 2:ny+1))
    
        for v in 1:k_iters
            # compute norm of laplacian error
            # lap_p_err_v = @view ARRAYS.lap_p_err[v,:, :]
            # err_view    = @view ARRAYS.err[2:nx+1, 2:ny+1]
            
            tmp_den[v] = sum(ARRAYS.lap_p_err[v,:, :].^2)
            if tmp_den[v] != 0.0
                beta = - dot(vec(ARRAYS.lap_p_err[v,:, :]), vec(ARRAYS.err[idx_interior])) / tmp_den[v]
            else
                beta = 0.25
            end
    
            # update Phi and error
            ARRAYS.Phi[idx_interior] .+= beta .* ARRAYS.p_err[v,idx_interior]
            ARRAYS.err[idx_interior] .+= beta .* ARRAYS.lap_p_err[v,:,:]

            xchng_pres_old!(ARRAYS.Phi)
    
            # check convergence
            error_step = maximum(abs.(ARRAYS.err[idx_interior]))
            if error_step <= error_tol
                converged = true
            end
    
            # compute new laplacian of error
            xchng_pres_old!(ARRAYS.err)
            ARRAYS.lap_err .= laplacian(ARRAYS.err, false, ARRAYS, dx, dy)
    
            # compute alpha coefficients
            for l in 1:v-1
                if tmp_den[l] != 0.0
                    # lap_p_err_l = @view ARRAYS.lap_p_err[l,:,:]
                    alpha[l] = - dot(vec(ARRAYS.lap_err),vec(ARRAYS.lap_p_err[l,:,:])) / tmp_den[l]
                end
            end
    
            # update p_err and lap_p_err for next iteration
            if v < k_iters
                ARRAYS.p_err[v+1,idx_interior] .= ARRAYS.err[idx_interior]
                ARRAYS.lap_p_err[v+1,:,:] .= ARRAYS.lap_err
                for l in 1:v-1
                    ARRAYS.p_err[v+1,idx_interior] .+= alpha[l] .* ARRAYS.p_err[l,idx_interior]
                    ARRAYS.lap_p_err[v+1,:,:] .+= alpha[l] .* ARRAYS.lap_p_err[l,:,:]
                end
            else
                # wrap around to first index
                ARRAYS.p_err[1,idx_interior] .= ARRAYS.err[idx_interior] .+ alpha[1] .* ARRAYS.p_err[1,idx_interior]
                ARRAYS.lap_p_err[1,:,:] .= ARRAYS.lap_err .+ alpha[1] .* ARRAYS.lap_p_err[1,:,:]
                for l in 2:v
                    ARRAYS.p_err[1,idx_interior] .+= alpha[l] .* ARRAYS.p_err[l,idx_interior]
                    ARRAYS.lap_p_err[1,:,:] .+= alpha[l] .* ARRAYS.lap_p_err[l,:,:]
                end
            end
        end
    
        return converged
end


function xchng_pres_old!(matrix::AbstractMatrix, halo::Int=1)
    matrix[1:halo, :] .= matrix[end-2*halo+1:end-halo, :]
    matrix[end-halo+1:end, :] .= matrix[halo+1:2*halo, :]

    matrix[:, 1:halo] .= matrix[:, end-2*halo+1:end-halo]
    matrix[:, end-halo+1:end] .= matrix[:, halo+1:2*halo]
end

function pressure_loop_init(ARRAYS, dx, dy)
    idx_interior = CartesianIndices((2:nx+1, 2:ny+1))
    ARRAYS.p_err[1,idx_interior] .= ARRAYS.err[idx_interior]
    xchng_pres_old!(@view ARRAYS.p_err[1, :,:])
    ARRAYS.lap_p_err[1,:,:] .= laplacian(ARRAYS.p_err[1,:,:], false,ARRAYS, dx, dy)
end

function pressure_solver_update(ARRAYS,u,v)


    idx_interior = CartesianIndices((2:nx+1, 2:ny+1))
    ARRAYS.tmp_uvw[1,idx_interior] .= u 
    ARRAYS.tmp_uvw[2,idx_interior] .= v

    xchng_pres_old!(@view ARRAYS.tmp_uvw[1, :,:]) 
    xchng_pres_old!(@view ARRAYS.tmp_uvw[2, :,:])
    ARRAYS.err[idx_interior] .= laplacian(ARRAYS.Phi, true,ARRAYS, dx, dy)
    xchng_pres_old!(ARRAYS.err)
    pressure_loop_init(ARRAYS, dx, dy)

    iters, converged = 0, false
    while converged == false      
        converged = pressure_solver_loop_body!(ARRAYS,converged)
        iters += 1
        if iters > 10000
            error(":(")
        end
    end
    xchng_pres_old!(ARRAYS.Phi)
    calc_grad_phi(ARRAYS.tmp_uvw, ARRAYS.Phi,idx_interior, dx, dy)
    xchng_pres_old!(@view ARRAYS.tmp_uvw[1,:,:])
    xchng_pres_old!(@view ARRAYS.tmp_uvw[2,:,:])

end

function pressure_solver_apply(ARRAYS,u,v,dt,nx,ny)
    idx_interior = CartesianIndices((2:nx+1, 2:ny+1))

    #does u and v really need to have halos?

    u .-= ARRAYS.tmp_uvw[1,idx_interior]
    ARRAYS.vip_rhs[1,idx_interior] .+= u
    ARRAYS.vip_rhs[1,idx_interior] ./= 0.5 * dt

    v .-= ARRAYS.tmp_uvw[2,idx_interior] 
    ARRAYS.vip_rhs[2,idx_interior] .+= v
    ARRAYS.vip_rhs[2,idx_interior] ./= 0.5 * dt

    xchng_pres_old!(@view ARRAYS.vip_rhs[1,:,:])
    xchng_pres_old!(@view ARRAYS.vip_rhs[2,:,:])
end

function vip_rhs_apply(ARRAYS, u,v,dt,nx,ny)
    idx_interior = CartesianIndices((2:nx+1, 2:ny+1))
    u .+= 0.5 * dt * ARRAYS.vip_rhs[1,idx_interior]
    ARRAYS.vip_rhs[1,:,:] .= 0
    v .+= 0.5 * dt * ARRAYS.vip_rhs[2,idx_interior]
    ARRAYS.vip_rhs[2,:,:] .= 0
end

function apply_half_rhs(ARRAYS,u,v,dt)              
    v .+= ARRAYS.rhs_w * dt/2
end 

function fill_stash(ARRAYS,u,v)
    idx_interior = CartesianIndices((2:nx+1, 2:ny+1))
    ARRAYS.stash[1,idx_interior] .= u
    ARRAYS.stash[2,idx_interior] .= v
    xchng_pres_old!(@view ARRAYS.stash[1, :,:])
    xchng_pres_old!(@view ARRAYS.stash[2, :,:])
end

function calc_gc_extrapolate_in_time(ARRAYS, u,v)
    idx_interior = CartesianIndices((2:nx+1, 2:ny+1))
    for _ in 1:2
        ARRAYS.stash[1,idx_interior] .= -.5 .* ARRAYS.stash[1,idx_interior] .+ 1.5 .* u
        ARRAYS.stash[2,idx_interior] .= -.5 * ARRAYS.stash[2,idx_interior] .+ 1.5 .* v
        xchng_pres_old!(@view ARRAYS.stash[1,:,:])
        xchng_pres_old!(@view ARRAYS.stash[2,:,:])
    end
end


function calc_gc_interpolate_in_space(ARRAYS, GCx, GCy, dx,dy,dt) 

    xchng_pres_old!(@view ARRAYS.stash[1,:,:])
    xchng_pres_old!(@view ARRAYS.stash[2,:,:])

    GCx .=  dt/dx *(
            (
                ARRAYS.stash[1,2:end,2:end-1] + 
                ARRAYS.stash[1,1:end-1,2:end-1]
            ) / 2
        )

    GCy .= dt/dy *(
            (
                ARRAYS.stash[2,2:end-1,2:end] + 
                ARRAYS.stash[2,2:end-1,1:end-1]
            ) / 2
        )

    #make sure first and last rows/columns are periodic
    GCx[1,:] .= GCx[end,:]
    GCx[end,:] .= GCx[1,:]
    GCy[:,1] .= GCy[:,end]
    GCy[:,end] .= GCy[:,1]
end

function finilize_rhs(ARRAYS,u,v,nx,ny)
    idx_interior = CartesianIndices((2:nx+1, 2:ny+1))
    ARRAYS.vip_rhs[1,idx_interior] .= -u
    ARRAYS.vip_rhs[2,idx_interior] .= -v
    xchng_pres_old!(@view ARRAYS.vip_rhs[1,:,:])
    xchng_pres_old!(@view ARRAYS.vip_rhs[2,:,:])
end

struct pressure_loop_tmp_struct
    rhs_w::Array{Float64,2}
    stash::Array{Float64,3}
    vip_rhs::Array{Float64,3}
    Phi::Array{Float64,2}
    tmp_uvw::Array{Float64,3}
    lap_tmp::Array{Float64,3}
    lap_err::Array{Float64,2}
    err::Array{Float64,2}
    p_err::Array{Float64,3}
    lap_p_err::Array{Float64,3}

    function pressure_loop_tmp_struct()
        rhs_w = zeros((nx, ny))
        stash = zeros((2, nx+2, ny+2))
        vip_rhs = zeros((2,nx+2, ny+2))
        Phi = zeros((nx+2, ny+2))
        tmp_uvw = zeros((2, nx+2, ny+2))
        lap_tmp = zeros((2,nx+2, ny+2))
        lap_err = zeros((nx, ny))
        err = zeros((nx+2, ny+2))
        p_err = zeros((4,nx+2, ny+2))
        lap_p_err = zeros((4,nx, ny))
        new(rhs_w,stash,vip_rhs,Phi,tmp_uvw,lap_tmp,lap_err,err,p_err,lap_p_err)
    end
end



