#---------------------------------------------------------
# SAMPLING
#---------------------------------------------------------
export init_ξ_const,init_logarithmic,init_uniform_sd, init_monodisperse, init_spatial_uniformvol_dry_sd

"""
    init_ξ_const(settings::coag_settings{FT}) where FT<:AbstractFloat

init_ξ_const initializes droplets based on the constant-multiplicity initialization method,
using an exponential distribution around the initial volume of the droplets. The method is
as described in Shima et al. (2009) https://doi.org/10.5194/acp-9-4491-2009

Arguments
- `settings`: Coagulation settings, type ::coag_settings.

"""
function init_ξ_const(settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    ξstart::Vector{Int} = (div(n0*ΔV,Ns)*ones(Ns))
    X0 = radius_to_volume(R0)# initial volume m3    
    Xstart::Vector{FT} = (rand(Exponential(X0), Ns))

    droplets = droplet_attributes(ξstart, Xstart)
    return droplets
end

"""
    init_logarithmic(settings::coag_settings{FT}) where FT<:AbstractFloat

init_logarithmic initializes the superdroplets, using a logarithmically spaced droplet radius spectrum,
initializing the multiplicities so that the water volume in the system forms an exponential distribution
around the initial volume.

Arguments
- `settings`: Coagulation settings.

"""

function init_logarithmic(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    R_min = settings.R_min
    R_max = settings.R_max

    X0 = radius_to_volume(R0)
    exp_dist = Exponential(X0)
    boundaries_found = false
    while boundaries_found == false
        radius_bins = 10 .^ range(log10(R_min), log10(R_max), length=Ns+1)
        # Calculate the volume for each bin
        vmin = radius_to_volume(radius_bins[1])
        vmax = radius_to_volume(radius_bins[end-1])
        
        pdf_min = pdf(exp_dist, vmin)
        pdf_max = pdf(exp_dist, vmax)

        first_bin_dr = radius_bins[2]-radius_bins[1]
        last_bin_dr = radius_bins[end]-radius_bins[end-1]
        dvdr_first = 4 * π * radius_bins[1].^2
        dvdr_last = 4 * π * radius_bins[end-1].^2

        ξ_first = pdf_min * first_bin_dr * dvdr_first * (n0*ΔV)
        ξ_last = pdf_max * last_bin_dr * dvdr_last * (n0*ΔV)

        if ξ_first >= 1 && ξ_last >=1
            boundaries_found = true
        else 
            if ξ_first <= 1
                R_min *= 1.01
            end
            if ξ_last <= 1
                R_max *= 0.99
            end
        end
    end

    radius_bins_new = 10 .^ range(log10(R_min), log10(R_max), length=Ns+1)
    sd_radii::Vector{FT} = [rand(Uniform(radius_bins_new[i], radius_bins_new[i+1])) for i in 1:Ns]
    volumes::Vector{FT} = radius_to_volume.(sd_radii)
    pdf_values = pdf.(exp_dist, volumes)
    bin_widths_new = diff(radius_bins_new)
    dvdr = 4 * π .* sd_radii.^2
    multiplicities = pdf_values .* bin_widths_new .* dvdr * (n0*ΔV)
    ξstart::Vector{Int} = floor.(multiplicities.+0.5)

    droplets = droplet_attributes(ξstart, volumes)
    return droplets
end

"""
init_uniform_sd(settings::coag_settings{FT}) where FT<:AbstractFloat

init_uniform_sd initializes the superdroplets, using an evenly spaced droplet radius spectrum,
initializing the multiplicities so that the water volume in the system forms an exponential distribution
around the initial volume.

Arguments
- `settings`: Coagulation settings.

"""

function init_uniform_sd(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    R_min = settings.R_min
    R_max = settings.R_max

    X0 = radius_to_volume(R0)
    exp_dist = Exponential(X0)
    boundaries_found = false
    while boundaries_found == false
        radius_bins = range(R_min, R_max, length=Ns+1)
        # Calculate the volume for each bin
        vmin = radius_to_volume.(radius_bins[1])
        vmax = radius_to_volume.(radius_bins[end-1])
        
        pdf_min = pdf(exp_dist, vmin)
        pdf_max = pdf(exp_dist, vmax)

        first_bin_dr = radius_bins[2]-radius_bins[1]
        last_bin_dr = radius_bins[end]-radius_bins[end-1]
        dvdr_first = 4 * π * radius_bins[1].^2
        dvdr_last = 4 * π * radius_bins[end-1].^2

        ξ_first = pdf_min * first_bin_dr * dvdr_first * (n0*ΔV)
        ξ_last = pdf_max * last_bin_dr * dvdr_last * (n0*ΔV)

        if ξ_first >= 1 && ξ_last >=1
            boundaries_found = true
        else 
            if ξ_first <= 1
                R_min *= 1.01
            end
            if ξ_last <= 1
                R_max *= 0.99
            end
        end
    end

    radius_bins_new = range(R_min, R_max, length=Ns+1)
    sd_radii::Vector{FT} = [rand(Uniform(radius_bins_new[i], radius_bins_new[i+1])) for i in 1:Ns]
    volumes::Vector{FT} = radius_to_volume.(sd_radii)
    pdf_values = pdf.(exp_dist, volumes)
    bin_widths_new = diff(radius_bins_new)
    dvdr = 4 * π .* sd_radii.^2
    multiplicities = pdf_values .* bin_widths_new .* dvdr * (n0*ΔV)
    ξstart::Vector{Int} = floor.(multiplicities.+0.5)

    droplets = droplet_attributes(ξstart, volumes)
    return droplets
end

"""
    init_monodisperse(settings::coag_settings{FT}) where FT<:AbstractFloat

init_monodisperse initializes the superdroplets so that all droplets have the same attributes.
Arguments
- `settings`: Coagulation settings.

"""

function init_monodisperse(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    ξstart::Vector{Int} = (div(n0*ΔV,Ns)*ones(Ns))
    Xstart::Vector{FT} = radius_to_volume(R0).*ones(Ns)

    droplets = droplet_attributes(ξstart,Xstart)
    return droplets
end


function init_spatial_uniformvol_dry_sd(rad_dist,settings::coag_settings{FT},
    spatial::spatial_settings{FT})where FT<:AbstractFloat

    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0

    percentile_limit = [0.001, 0.999]  
    Rmin = cdf(rad_dist, percentile_limit[1])
    Rmax = cdf(rad_dist, percentile_limit[2])   

    Rarray = range(Rmin, Rmax, length=Ns+1)
    cdf_values = cdf.(dist, Rarray)
    rad = (Rarray[2:end] .+ Rarray[1:end-1])./ 2 
    cdf_values = cdf_values[2:end] - cdf_values[1:end-1]
    multiplicities = cdf_values .* (n0 * ΔV)
                    
    ξstart::Vector{Int} = floor.(multiplicities .+ 0.5)
    z_loc_in_cell = rand(Uniform(0, 1), Ns)
    x_loc_in_cell = rand(Uniform(0, 1), Ns)
    cell_id = rand(1:spatial.num_grids, Ns)

    dry_r3 = rad.^3

    #set to activation radius
    r_crit = sqrt(3 * B / A)/ 2
    volumes = 4* pi / 3 .*dry_r3

    droplets = droplet_attributes_2d{FT}(ξstart, volumes,dry_r3,z_loc_in_cell, x_loc_in_cell, cell_id)
    return droplets
end
