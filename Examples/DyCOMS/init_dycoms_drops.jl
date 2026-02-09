using Distributions
using Droplets
using Plots

function init_droplets_dycoms_scm(settings::coag_settings{FT},
    spatial::spatial_settings{FT})where FT<:AbstractFloat

    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0

    #double The total number,
    # mode radius, and geometric standard deviation for the
    # two modes are 125 and 65 cm23, 0.011 and 0.06 um, and
    # 1.2 and 1.7,
    m1 = 0.011e-6 # mode 1 in meters
    σ1 = log(1.2) # geometric standard deviation for mode 1
    m2 = 0.06e-6 # mode 2 in meters
    σ2 = log(1.7) # geometric standard deviation for mode 2
    μ1 = log(m1) + σ1^2 
    μ2 = log(m2) + σ2^2
    p1 = 125 / (125 + 65) 
    p2 = 65 / (125 + 65) 
    dist = MixtureModel(LogNormal, [(μ1,σ1),(μ2,σ2)], [p1, p2])



    percentile_limit = [0.001, 0.999]  
    Rmin = quantile(dist, percentile_limit[1])
    Rmax = quantile(dist, percentile_limit[2])   

    Rarray = range(Rmin, Rmax, length=Ns+1)
    cdf_values = cdf.(dist, Rarray)
    rad = (Rarray[2:end] .+ Rarray[1:end-1])./ 2 
    cdf_values = cdf_values[2:end] - cdf_values[1:end-1]
    multiplicities = cdf_values .* (n0 * spatial.Z_max * spatial.area_per_grid)
                    
    ξstart::Vector{Int} = floor.(multiplicities .+ 0.5)
    z_loc = rand(Uniform(0, 1), Ns) * spatial.Z_max
    cell_id = floor.(z_loc ./ spatial.z_grid_height) .+ 1

    dry_r3 = rad.^3
    kappa = 0.61 * ones(Ns) #ammonium sulfate, Petters and Kreidenweis (2007)

    volumes = 4* pi / 3 .* dry_r3 #initialize as just aerosol volume.

    droplets = droplet_attributes_1d{FT}(ξstart, volumes,dry_r3,kappa,z_loc, cell_id)
    return droplets
end