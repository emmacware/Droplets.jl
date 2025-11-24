using Plots

function coag_runtime(randseed::Int,droplets::droplet_attributes,
    coag_settings::coag_settings{FT},run_settings::run_settings{FT}) where FT<:AbstractFloat
    
    Random.seed!(randseed)
    println("Running simulation...")

    coal_func_time::FT = 0.0
    bins::Matrix{FT} = zeros(FT, run_settings.num_bins, length(run_settings.output_steps))
    threading,scheme = run_settings.coag_threading, run_settings.scheme
    coag_data = coagulation_run{FT}(coag_settings.Ns)
    simtime::FT = @CPUelapsed begin
        for i  in  1:length(run_settings.output_steps)
            # if i,seconds in enumerate(run_settings.output_steps)
            
            if i !=1
                timestepper = Int(round((run_settings.output_steps[i]-run_settings.output_steps[i-1])/coag_settings.Δt))
                ctime::FT = @CPUelapsed begin
                    for _ in 1:timestepper
                        coalescence_timestep!(threading,scheme,droplets,coag_data,coag_settings)
                    end
                end
                coal_func_time += ctime
            end
            bins[:,i] = binning_func(droplets,run_settings.output_steps[i],run_settings,coag_settings)
            # println("Time: ", run_settings.output_steps[i], " seconds")
        end
    end
    println("simtime =", simtime)
    println("coal_func_time =", coal_func_time)

    return bins, coal_func_time
end


function plot_dsd(bins,runsettings::run_settings{FT};color="black",label=false,legend=true) where FT<:AbstractFloat   
    radius_bins_edges = runsettings.radius_bins_edges
    mids = 0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])*1e6
    plot1 = plot!(mids,bins,xaxis=:log,label= runsettings.output_steps',legend=legend,legendtitle="Time Steps (s)")
    xlabel!("Radius (μm)")
    if runsettings.binning_method == number_density
        ylabel!("Number density (m^-3)")
    elseif runsettings.binning_method == mass_density_lnr
        ylabel!("Mass density (kg m^-3)")
    end
    title!("Time Evolution of Droplet Size Distribution")

    return plot1
end
