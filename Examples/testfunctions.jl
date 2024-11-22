using Plots



# function build_box_run(coag_settings::coag_settings{FT}, run_settings::run_settings{FT}) where FT<:AbstractFloat
#     ξ, R, X = run_settings.init_method(coag_settings)
#     Ns::Int = coag_settings.Ns
#     I::Vector{Int} = shuffle(1:Ns)
    
#     function coag_runtime(randseed::Int,ξ::Vector{Int},R::Vector{FT},X::Vector{FT},
#         I::Vector{Int},coag_settings::coag_settings{FT},run_settings::run_settings{FT}) where FT<:AbstractFloat
#         Random.seed!(randseed)
        
#         println("Running simulation...")
#         seconds::FT = 0.0

#         coal_func_time::FT = 0.0
#         bins::Matrix{FT} = zeros(FT, run_settings.num_bins - 1, 0)

#         ϕ = Vector{FT}(undef, div(coag_settings.Ns, 2))
#         simtime::FT = @elapsed begin
#             while seconds <= run_settings.output_steps[end]
#                 if seconds in run_settings.output_steps
#                     println("Time: ", seconds, " seconds")
#                     xx, yy = run_settings.binning_method(X, ξ,seconds, run_settings)
#                     bins = hcat(bins, yy)
#                 end

#                 ctime::FT = @elapsed begin
#                     coalescence_timestep!(run_settings.coag_threading, run_settings.scheme, ξ, R, X, I,ϕ,coag_settings)
#                 end
#                 coal_func_time += ctime
#                 seconds += coag_settings.Δt
#                 # GC.gc()
#             end
#         end
#         println("simtime =", simtime)
#         println("coal_func_time =", coal_func_time)

#         return bins, coal_func_time
#     end
#     return coag_runtime
# end

function coag_runtime(randseed::Int,droplets::droplets_allocations,
    coag_settings::coag_settings{FT},run_settings::run_settings{FT}) where FT<:AbstractFloat
    
    Random.seed!(randseed)
    println("Running simulation...")

    coal_func_time::FT = 0.0
    bins::Matrix{FT} = zeros(FT, run_settings.num_bins - 1, length(run_settings.output_steps))
    threading,scheme = run_settings.coag_threading, run_settings.scheme
    simtime::FT = @CPUelapsed begin
        for i  in  1:length(run_settings.output_steps)
            # if i,seconds in enumerate(run_settings.output_steps)
            
            if i ==1
                bins[:,i] = run_settings.binning_method(droplets.X, droplets.ξ,run_settings.output_steps[i],run_settings)
                println("Time: ", run_settings.output_steps[i], " seconds")
                continue
            end

            timestepper = (run_settings.output_steps[i]-run_settings.output_steps[i-1])/coag_settings.Δt
            ctime::FT = @CPUelapsed begin
                for _ in 1:timestepper
                    coalescence_timestep!(threading,scheme,droplets,coag_settings)
                end
            end
            coal_func_time += ctime
            bins[:,i] = mass_density_lnr(droplets.X, droplets.ξ,run_settings.output_steps[i],run_settings)
            println("Time: ", run_settings.output_steps[i], " seconds")
        end
    end
    println("simtime =", simtime)
    println("coal_func_time =", coal_func_time)

    return bins, coal_func_time
end

