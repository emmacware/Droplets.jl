#---------------------------------------------------------
#Visualization functions for droplet size distribution
#---------------------------------------------------------
export run_settings, error_measure, number_density, mass_density_lnr, binning_func, binning_1d, smoothbins!

#####

"""
    struct run_settings{FT<:AbstractFloat}

A struct representing the run settings for the binning process.

Fields:
- `FT`: The type of floating-point numbers used in the run settings.
- `num_bins`: The number of bins to use for binning output(default: 128).
- `radius_bins_edges`: The edges of the bins for the radius (default: logarithmic scale from 10e-6 to 5e3e-6).
- `smooth`: A boolean indicating whether to apply smoothing to the bins (default: true).
- `smooth_scope`: The scope of smoothing (default: 2).
- `init_random_seed`: The random seed for initialization (default: 30).
- `coag_threading`: The threading method for coalescence (default: Serial).
    - Serial() or Parallel()
- `scheme`: The scheme to use for coalescence (default: none).
    -none(), Adaptive()
- `output_steps`: The time steps at which to output the results (default: [0, 1200, 2400, 3600]).
- `init_method`: The method to use for initialization (default: init_logarithmic).
    - init_logarithmic,init_ξ_const, init_uniform_sd
- `binning_method`: The method to use for binning (default: mass_density_lnr). 
    (a function that takes droplets and coagsettings as arguments)
- `normalize_bins_dlnr`: A boolean indicating whether to normalize the bins by the logarithm of the radius (default: true).

"""
Base.@kwdef struct run_settings{FT<:AbstractFloat} #??
    num_bins::Int = Int(128)
    radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=num_bins+1) 
    smooth::Bool = true #bin smoothing
    smooth_scope::Int = Int(2)
    init_random_seed::Int = Int(30) 
    coag_threading =  Serial()#Parallel(),use Julia NThreads for coalescence
    scheme = none() #Adaptive,Small_Alpha
    output_steps::Vector{FT} = [0,1200,2400,3600]
    init_method = init_logarithmic #init_ξ_const,init_uniform_sd
    # spectrum = exponential #only this for now
    binning_method = mass_density_lnr #number_density
    normalize_bins_dlnr = true
    # condensation::string = "none" #not yet
    #decide output here? bins,error,timing,...
end


#integrated error function, as in PySDM
function error_measure_old(y,ytrue,x)
    errors = ytrue .- y
    errors1 = errors[1:end-1] .+ errors[2:end]
    dx = diff(x)
    errors1 = errors1 .* dx ./2
    err = sum(abs.(errors1))
    return err
end

function error_measure(y,ytrue)
    return sqrt(sum((y.-ytrue).^2)/length(y))
end


function number_density(droplets::droplet_attributes,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    weights = droplets.ξ/coagsettings.ΔV
    return weights
end


function mass_density_lnr(droplets::droplet_attributes,coagsettings::coag_settings{FT}) where FT<:AbstractFloat

    tot_vol = droplets.ξ .* droplets.X
    weights_kilograms = tot_vol*constants.ρl # convert to mass
    kg_per_vol = weights_kilograms/coagsettings.ΔV 

    return kg_per_vol
end

function binning_func(droplets::droplet_attributes, t::FT,
    runsettings::run_settings{FT}, coagsettings::coag_settings{FT};
    indices = eachindex(droplets.X)) where FT<:AbstractFloat

    weights = runsettings.binning_method(droplets,coagsettings)
    R = volume_to_radius.(droplets.X)
    idx = droplets.I[indices]
    numdens = binning_1d(R[idx], weights[idx], runsettings)

    if t != 0 && runsettings.smooth == true
        numdens = smoothbins!(numdens,runsettings)
    end

    return numdens
end

function binning_1d(r_values_unsorted::Vector{},weights_unsorted::Vector{},runsettings::run_settings{FT}) where FT<:AbstractFloat
    bin_edges = runsettings.radius_bins_edges

    numdens::Vector{FT} = zeros(runsettings.num_bins)
    idx = sortperm(r_values_unsorted)
    values = r_values_unsorted[idx]
    weights = weights_unsorted[idx]

    droplet_idx = 1
    for j in 1:runsettings.num_bins
        bin_edge_high = (bin_edges[j+1])
        for i in droplet_idx:length(values)
            if values[i] < bin_edge_high
                numdens[j]=numdens[j]+weights[i]
                droplet_idx += 1
            else
                break
            end
        end
    end

    if runsettings.normalize_bins_dlnr == true
        numdens = numdens ./ diff(log.(bin_edges))
    else
        numdens = numdens ./ diff(bin_edges)
    end

    return numdens
end



function smoothbins!(numdens::Vector{FT},runsettings::run_settings{FT}) where FT<:AbstractFloat
    scope = runsettings.smooth_scope

    new_numdens = copy(numdens)
    for _ in 1:2

        for i in scope+1:length(numdens)-scope
            new_numdens[i] = sum(numdens[i-scope:i+scope])/length(numdens[i-scope:i+scope])
        end
        scope = 1

        for i in scope+1:length(numdens)-scope
            numdens[i] = sum(new_numdens[i-scope:i+scope])/length(new_numdens[i-scope:i+scope])
        end
    end
    return numdens
end


