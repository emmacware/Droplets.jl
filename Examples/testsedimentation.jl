using Random
using Distributions
using Plots
using Interpolations
using DifferentialEquations
include("DSDvis.jl")
using Droplets



# ##################################################
# Sedimentation example, combining condensation and coalescence,
# allowing for vertical movement of droplets

# ##################################################

# #Settings 

#spatial domain (single column)
#-------------------------------------------------
Nx = 1
Ny = 5
Δx = 100
Δy = 100
qv = [0.01, 0.01, 0.01, 0.01, 0.01] #kg/kg
P = [101.3, 101.3, 101.3, 101.3, 101.3] #hPa

#-------------------------------------------------
# Grid
grid_box, grid_box_mids_x, grid_box_mids_y = create_gridbox(Nx,Ny,Δx,Δy)

# #-------------------------------------------------
Sv = zeros(Nx,Ny)
P0 = 101325 #Pa
T = zeros(Nx,Ny).+305
for j in 1:Ny
    height = grid_box_mids_y[1,j]/1000
    T[:,j] .= 305-9.5*height
end

# sed_const = Constants()

qvarray = [0.0295 0.028 0.027 0.026 0.025] #kg/kg, specific humidity or mixing ratio of vapor/moist air
P = P0.*exp.(-constants.gconst.*grid_box_mids_y./(T.*constants.Rd))
ρu = zeros(Nx,Ny) #just sedimentation
ρv = zeros(Nx,Ny+1) #just sedimentation
ρd =  P./(constants.Rd.*T)
ρ = ρd ./(1 .- qvarray) #kg/m^3, density of moist air
θb = T.*(P0./P).^(constants.Rd/constants.Cp)

#superdroplets
# #-------------------------------------------------
Ns = 2^8
n0 = 2^23
R0 = 30.531e-6 #meters
M0 = 1e-16 #kg
Δt = 1
ΔV = 10^2*Nx*Δx*Ny*Δy #depth of 100m, full system
kernel = golovin
radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128) 
smooth = true

ξstart,Rstart,Xstart,Mstart = init_ξ_const(Ns,ΔV,n0,R0,M0)

superdroplets = create_NaCl_superdroplets(Ns,Nx,Ny,Δx,Δy,Rstart,Xstart,ξstart,Mstart)
grid_dict = Dict{Tuple{Int, Int}, Vector{Superdroplet}}()
droplet_gridbox(superdroplets,Nx,Ny,Δx,Δy,grid_dict)

# #-------------------------------------------------
function sedimentation_animation(superdroplets,grid_dict,grid_box,grid_box_mids_y,Δt,ΔV,Nx,Ny,ρu,ρv,ρ,qvarray,P,T,θb,radius_bins_edges,kernel)
anim = Animation()
plot3 = heatmap([1],grid_box_mids_y[:],T',c=cgrad(:roma, rev = true),colorbar_title="T(K)")
Senv = sat.(qvarray,P)./esat.(T)
plot3 = plot!(Senv',grid_box_mids_y[:],c=:black,linewidth=4,label = "Senv",background=:dodgerblue)
xlims!(plot3,0.995,1.025)
xlabel!(plot3,"Saturation")
title!(plot3,"T,Senv")
xticks!(plot3,[1.0,1.01,1.02],["1.0","1.01","1.02"])
plot2 = plot()
plot2 = plot_grid_with_droplets(grid_box, superdroplets)
plot2 = ylims!(0,500)
plot2 = xlabel!("X (m)")
plot2 = ylabel!("Z (m)")
plot2 = title!("Δt = 1m,Sedimentation, Coalescence, and Condensation",titlefontsize=12)
# title1 = plot(title = "Sedimentation, Coalescence, and Condensation", grid = false, showaxis = false)

plot1 = [plot(),plot(),plot(),plot(),plot()]
for i in 1:Ny
            droplets = grid_dict[(1,i)]
            ngrid = length(droplets)
            x,y = binning_dsd(droplets,radius_bins_edges,Δt;smooth = true,scope_init = 2)
            # plot1[i]= plot(x,y,log)
            plot1[Ny+1-i] = plot(x,y,lc=:black,xaxis=:log,xlims=(10,1000),ylims=(0,3),label=false)
end
xlabel!(plot1[Ny], "Radius (μm)")
ylabel!(plot1[3], "g(lnr)")
title!(plot1[1],"DSD")
p1 = plot(plot1...,layout=(Ny,1))
plot(plot3,plot2,p1,layout=@layout([A B{0.5w} C]),size=(1000,600))
# plot(plot2, plot1..., layout=@layout([a{0.7w} grid(Ny, 1)]))
# ylims!(0,500)
frame(anim)
sleep(1)

for t in 1:15
    for _ in 1:60
        for i in 1:Ny
            droplets = grid_dict[(1,i)]
            ngrid = length(droplets)
            coalescence_timestep!(droplets,ngrid,Δt,ΔV,kernel=kernel)
        #     x,y = binning_dsd(droplets,radius_bins_edges,Δt;smooth = true,scope_init = 2)
        #     # plot1[i]= plot(x,y,log)
        #     plot1[Ny+1-i] = plot(x,y,lc=:black,xaxis=:log,xlims=(10,1000),ylims=(0,3),label=false)
        end
        Sv = condense_and_calc_Sv!(qvarray,T,P,ρ,Δt,ΔV,Nx,Ny,grid_dict)
        θb,T = θcondenseupdate!(Sv,θb,Δt,P,P0)
        qvarray,ρ = qvcondenseupdate!(Sv, qvarray, P,T, Δt)
        update_position!(superdroplets,Nx,Ny,Δt,ρu,ρv,ρ, grid_dict,grid_box)
    end
    for i in 1:Ny
        droplets = grid_dict[(1,i)]
        ngrid = length(droplets)
        x,y = binning_dsd(droplets,radius_bins_edges,Δt;smooth = true,scope_init = 2)
        # plot1[i]= plot(x,y,log)
        plot1[Ny+1-i] = plot(x,y,lc=:black,xaxis=:log,xlims=(10,1000),ylims=(0,3),label=false)
    end
    plot3 = heatmap([1],grid_box_mids_y[:],T',c=cgrad(:roma, rev = true),colorbar_title="T(K)")
    Senv = sat.(qvarray,P)./esat.(T)
    plot3 = plot!(Senv',grid_box_mids_y[:],c=:black,linewidth=4,label = "Senv",background=:dodgerblue)
    title!(plot3,"T,Senv")
    xlims!(plot3,0.995,1.025)
    xlabel!(plot3,"Saturation")
    xticks!(plot3,[1.0,1.01,1.02],["1.0","1.01","1.02"])    
    plot2 = plot()
    plot2 = plot_grid_with_droplets(grid_box, superdroplets)
    plot2 = ylims!(0,500)
    plot2 = xlabel!("X (m)")
    plot2 = ylabel!("Z (m)")
    plot2 = title!("Δt = 1m,Sedimentation, Coalescence, and Condensation",titlefontsize=12)
    title!(plot1[1],"DSD")
    xlabel!(plot1[Ny], "Radius (μm)")
    ylabel!(plot1[3], "g(lnr)")
    p1 = plot(plot1...,layout=(Ny,1))
    plot(plot3,plot2,p1,layout=@layout([A B{0.5w} C]),size=(1000,600))

    frame(anim)
    sleep(1)
end
return anim
end

anim = sedimentation_animation(superdroplets,grid_dict,grid_box,grid_box_mids_y,Δt,ΔV,Nx,Ny,ρu,ρv,ρ,qvarray,P,T,θb,radius_bins_edges,kernel)
gif(anim, fps = 2)