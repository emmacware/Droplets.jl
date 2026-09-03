import CairoMakie as MK
import CloudMicrophysics as CM
using Droplets
using ComponentArrays

# Runs CM Liquid-only parcel example for initial conditions/variables and plot objects loaded
include(joinpath(pkgdir(CM), "parcel", "Example_Liquid_only.jl"))

# Add Droplets.jl <-> CM coupling (SDCondensation, parcel_model_sd, run_parcel_sd)
include("CM_sdfunctions.jl")

Ns = Int(2^14) #number of superdroplets

#Droplets has Exponential and Monodisperse, replace Gamma with Exponential for the SD runs
liq_size_distribution_list[findfirst(isequal("Gamma"), liq_size_distribution_list)] = "Exponential"


for DSD in liq_size_distribution_list
    local pp = parcel_params{FT}(
        liq_size_distribution = DSD,
        condensation_growth = "Superdroplets",
        const_dt = const_dt,
        w = w,
    )
    coagsettings = coag_settings{FT}(Ns = Ns, Δt = FT(1), ΔV = FT(1), n0 = Nₗ, R0 = r₀)
    local drops = DSD == "Monodisperse" ? init_monodisperse(coagsettings) : init_ξ_const(coagsettings)

    local ml_v = sum(drops.X .* drops.ξ .* ρₗ)
    local qᵥ = mv_v / (md_v + mv_v + ml_v)
    local qₗ = ml_v / (md_v + mv_v + ml_v)
    local IC_sd = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]
    local d_Y = ComponentVector{FT}(IC = IC_sd, ξ = drops.ξ, X = drops.X)

    local sol = run_parcel_sd(d_Y, FT(0), t_max, pp)

    label = DSD * " (SD)"
    MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0, label = label, linestyle = :dash)
    MK.lines!(ax3, sol.t, sol[5, :] * 1e3, linestyle = :dash)

    local r = similar(sol[3, :])
    for (i, output) in enumerate(sol.u)
        r[i] = (3 / (4π) * sum(output.X) / Ns)^(1 / 3)
    end
    MK.lines!(ax2, sol.t, r * 1e6, linestyle = :dash)
end

# # The included example already drew a legend with just the CLiMA curves -- drop it and
# # redraw with all four now that the SD curves are on the same axes.
for leg in filter(x -> x isa MK.Legend, fig.content)
    MK.delete!(leg)
end
MK.axislegend(ax1, framevisible = false, labelsize = 12, orientation = :vertical, position = :rb)

MK.save("liquid_only_parcel_coupled.svg", fig)
# nothing
