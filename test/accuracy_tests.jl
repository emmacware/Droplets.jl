using Test
using Droplets
using Random


FT = Float64

@testset "golovinSCE / density_lnr_golovin_analytic: SDM converges to the analytic Golovin solution" begin
    # compare w analytical solution
    Random.seed!(20240607)

    num_bins = 32
    radius_bins_edges = 10 .^ range(log10(1e-6), log10(5e-3), length=num_bins + 1)
    output_steps = FT[0, 1800, 3600]

    coagsettings = coag_settings{FT}(Ns=2^13, Δt=FT(1.0), ΔV=FT(1e6),
        golovin_kernel_coeff=FT(1500.0), n0=FT(2^23), R0=FT(30.531e-6), kernel=golovin)
    runsettings = run_settings{FT}(num_bins=num_bins, radius_bins_edges=radius_bins_edges,
        smooth=true, output_steps=output_steps, init_method=init_logarithmic,
        binning_method=mass_density_lnr)

    analytic = density_lnr_golovin_analytic(runsettings, coagsettings)
    @test all(isfinite, analytic)
    @test all(analytic .>= 0)

    droplets = runsettings.init_method(coagsettings)
    mass0 = sum(droplets.ξ .* droplets.X)
    coag_data = coagulation_run{FT}(coagsettings.Ns)

    sim_bins = zeros(FT, num_bins, length(output_steps))
    sim_bins[:, 1] = binning_func(droplets, output_steps[1], runsettings, coagsettings)

    for i in 2:length(output_steps)
        nsteps = Int(round((output_steps[i] - output_steps[i-1]) / coagsettings.Δt))
        for _ in 1:nsteps
            coalescence_timestep!(runsettings.coag_threading, runsettings.scheme, droplets, coag_data, coagsettings)
        end
        sim_bins[:, i] = binning_func(droplets, output_steps[i], runsettings, coagsettings)
    end

    # mass conservation throughout 
    @test sum(droplets.ξ .* droplets.X) ≈ mass0 rtol=1e-8

    # size-resolved comparison against the analytic Smoluchowski solution. SDM with a
    # finite (2^13) superdroplet count carries real Monte Carlo noise, so the tolerance
    # is generous
    for i in 2:length(output_steps)
        peak = maximum(analytic[:, i])
        err = error_measure(sim_bins[:, i], analytic[:, i])
        @test err < 0.5 * peak

        # the simulated and analytic spectra should peak in roughly the same size range
        analytic_peak_r = radius_bins_edges[argmax(analytic[:, i])]
        sim_peak_r = radius_bins_edges[argmax(sim_bins[:, i])]
        @test analytic_peak_r / sim_peak_r < 4 && sim_peak_r / analytic_peak_r < 4
    end

    # total droplet number strictly decreases under coalescence 
    @test sum(droplets.ξ) < coagsettings.n0 * coagsettings.ΔV
end

@testset "hydrodynamic kernel: symmetry, positivity, Hall-Davis efficiency wiring" begin
    coagsettings_noeff = coag_settings{FT}(hydrodynamic_collision_eff_func=false, kernel=hydrodynamic)
    coagsettings_eff = coag_settings{FT}(hydrodynamic_collision_eff_func=true, kernel=hydrodynamic)

    R = FT[1e-6, 5e-6, 20e-6, 100e-6, 500e-6]
    X = radius_to_volume.(R)
    drops = droplet_attributes{FT}(fill(1, length(R)), X)

    for (j, k) in [(1, 2), (1, 5), (3, 4), (2, 2)]
        for settings in (coagsettings_noeff, coagsettings_eff)
            Kjk = hydrodynamic(drops, (j, k), settings)
            Kkj = hydrodynamic(drops, (k, j), settings)
            @test Kjk ≈ Kkj  # kernel must be symmetric in its two arguments
            @test Kjk >= 0
            @test isfinite(Kjk)
        end
    end

    # equal-size droplets never collide (|v(r)-v(r)| = 0), efficiency notwithstanding
    @test hydrodynamic(drops, (2, 2), coagsettings_noeff) == 0.0
    @test hydrodynamic(drops, (2, 2), coagsettings_eff) == 0.0

    # turning on the efficiency function must actually change the kernel (E < 1
    # somewhere), or the flag is wired to nothing
    K_noeff = hydrodynamic(drops, (1, 3), coagsettings_noeff)
    K_eff = hydrodynamic(drops, (1, 3), coagsettings_eff)
    @test K_eff <= K_noeff
    K_noeff_big = hydrodynamic(drops, (3, 5), coagsettings_noeff)
    K_eff_big = hydrodynamic(drops, (3, 5), coagsettings_eff)
    @test K_eff_big <= K_noeff_big
    @test K_eff_big > 0
end

@testset "collision_efficiency: Hall-Davis table lookup sanity" begin
    # symmetry: swapping the two radii must give the same value,
    for (r1, r2) in [(5e-6, 30e-6), (12e-6, 12e-6), (60e-6, 250e-6), (600e-6, 900e-6)]
        @test Droplets.collision_efficiency(FT(r1), FT(r2)) ≈ Droplets.collision_efficiency(FT(r2), FT(r1)) atol=1e-10
    end

    # 1µm droplets: zero collision efficiency at the very bottom of the table
    @test Droplets.collision_efficiency(FT(1e-6), FT(1e-6)) == 0.0

    # all table lookups (interpolated or clamped) return finite, non-negative values
    for r1 in FT[1e-6, 3e-6, 15e-6, 50e-6, 150e-6, 600e-6, 1200e-6], r2 in FT[1e-6, 40e-6, 300e-6, 1500e-6]
        e = Droplets.collision_efficiency(r1, r2)
        @test isfinite(e)
        @test e >= 0
    end
end

@testset "number_density / mass_density_lnr" begin
    ξ = [10, 20, 30]
    R = FT[1e-6, 2e-6, 3e-6]
    X = radius_to_volume.(R)
    drops = droplet_attributes{FT}(ξ, X)
    coagsettings = coag_settings{FT}(ΔV=FT(2.0))

    nd = number_density(drops, coagsettings)
    @test nd ≈ ξ ./ 2.0

    md = mass_density_lnr(drops, coagsettings)
    @test md ≈ (ξ .* X .* constants.ρl) ./ 2.0
end

@testset "binning_1d: known radii land in the expected bins" begin
    bin_edges = FT[1.0, 2.0, 4.0, 8.0, 16.0]  # 4 bins
    runsettings = run_settings{FT}(num_bins=4, radius_bins_edges=bin_edges,
        smooth=false, normalize_bins_dlnr=false)

    r_values = FT[1.5, 1.9, 3.0, 3.5, 3.9, 10.0, 15.9]
    weights = ones(FT, length(r_values))

    numdens = binning_1d(r_values, weights, runsettings)

    @test numdens[1] ≈ 2.0 / (bin_edges[2] - bin_edges[1])   # 1.5, 1.9
    @test numdens[2] ≈ 3.0 / (bin_edges[3] - bin_edges[2])   # 3.0, 3.5, 3.9
    @test numdens[3] ≈ 0.0                                    # nothing in [4,8)
    @test numdens[4] ≈ 2.0 / (bin_edges[5] - bin_edges[4])   # 10.0, 15.9
    @test sum(numdens .* diff(bin_edges)) ≈ length(r_values)  # total weight conserved
end

@testset "binning_1d: dlnr normalization matches manual computation" begin
    bin_edges = FT[1.0, 10.0, 100.0]
    runsettings = run_settings{FT}(num_bins=2, radius_bins_edges=bin_edges,
        smooth=false, normalize_bins_dlnr=true)
    r_values = FT[5.0, 5.0, 50.0]
    weights = FT[1.0, 3.0, 2.0]

    numdens = binning_1d(r_values, weights, runsettings)
    @test numdens[1] ≈ 4.0 / log(10.0 / 1.0)
    @test numdens[2] ≈ 2.0 / log(100.0 / 10.0)
end

@testset "smoothbins!: leaves the two endpoint bins untouched" begin
    numdens = FT[1.0, 5.0, 2.0, 8.0, 3.0, 9.0, 1.0, 4.0, 2.0, 6.0]
    runsettings = run_settings{FT}(smooth_scope=2)
    original = copy(numdens)
    smoothbins!(numdens, runsettings)

    @test numdens[1] == original[1]
    @test numdens[end] == original[end]
    @test numdens[2:end-1] != original[2:end-1]  # interior bins actually get smoothed
    @test all(isfinite, numdens)
end
